#include "cosmology.h"
#include "ic2.h"
#include "interfaces.h"
#include "io.h"
#include "parameters.h"
#include "state.h"

/* #include <fftw3.h> */
#include <algorithm>
#include <fstream>
/* #include <map> */
#include <numeric>
/* #include <random> */
#include <vector>

using namespace blaze;

// init wavefunction ...
// by loading a file containing delta(x_i)
void delta_from_file(SimState&, std::ifstream&);
// init wavefunction from file
void psi_from_file(SimState&, std::ifstream&);
// according to a matter power spectrum provided by file.
void delta_from_power(SimState&, const Cosmology&, std::ifstream&);

void generate(SimState& state, const Cosmology& cosmo, const Parameters& p) {
    auto type =
        static_cast<ICType>(p["Initial Conditions"]["ic_type"].get<int>());
    std::ifstream ic_file{
        p["Initial Conditions"]["source_file"].get<std::string>()};

    switch (type) {
        case ICType::ExternalDelta... ICType::ExternalPsi: {
            size_t data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                                       std::istreambuf_iterator<char>(), '\n');
            if (data_N != state.N_total) {
                std::cout << ERRORTAG("#lines in source_file(" << data_N
                                                               << ") != N")
                          << std::endl;
                exit(1);
            }
            ic_file.seekg(0);
        }
        default:
            break;
    };

    switch (type) {
        case ICType::ExternalDelta:
            std::cout << INFOTAG("Load delta0 from file") << std::flush;
            delta_from_file(state, ic_file);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::ExternalPsi:
            std::cout << INFOTAG("Load Psi0 from file") << std::flush;
            psi_from_file(state, ic_file);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::Powerspectrum:
            std::cout << INFOTAG("Generate delta0 from P(k)") << std::flush;
            delta_from_power(state, cosmo, ic_file);
            std::cout << " ... done" << std::endl;
            break;
    }

    // Velocity Initialization
    if (type != ICType::ExternalPsi) {
        // state.V holds delta at this point. But we will use it in an in-place
        // manner. Hence we define...
        const auto& delta = state.V;
        auto& phase = state.V;
        auto& psi = state.psi;

        // Set psis modulus before we override delta. For cold conditions that's
        // all that is left to do.
        psi = sqrt(1.0 + delta);

        // Cosmological initial velocity field
        const bool compute_velocity =
            p["Initial Conditions"]["compute_velocity"].get<bool>();
        if (compute_velocity) {
            auto pot = Interaction::make("Poisson::FFT", p);
            const double a_init = cosmo.a_of_tau(0);
            const double prefactor =
                -std::sqrt(2.0 / 3 * a_init / cosmo.omega_m(a_init));

            std::cout << INFOTAG("Initial Velocity Field from Poisson @ a = ")
                      << a_init << std::endl;

            pot->solve(phase, prefactor * delta);

            // Madelung Representation.
            psi *= exp(std::complex<double>(0, 1) * phase);
        }
    }
}

void delta_from_file(SimState& state, std::ifstream& ic_file) {
    // Store delta in V for convenience (real vector vs complex vector);
    auto& delta = state.V;

    fill_from_file(ic_file, delta);
}

void psi_from_file(SimState& state, std::ifstream& ic_file) {
    DynamicVector<double> re(state.N_total);
    DynamicVector<double> im(state.N_total);

    fill_from_file(ic_file, re, im);
    state.psi =
        map(re, im, [](double r, double i) { return std::complex(r, i); });
}

// TODO Implement d+1 dimensional version of cosmological IC
void delta_from_power(SimState& state, const Cosmology& cosmo,
                      std::ifstream& ic_file) {
    exit(1);
    /*
    size_t data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                               std::istreambuf_iterator<char>(), '\n');
    ic_file.seekg(0);
    std::vector<double> k_data(data_N);
    std::vector<double> powerspectrum(data_N);
    std::map<double, double> discrete_P;
    fill_from_file(ic_file, k_data, powerspectrum);

    std::transform(k_data.begin(), k_data.end(), powerspectrum.begin(),
                   std::inserter(discrete_P, discrete_P.end()),
                   std::make_pair<const double&, const double&>);

    // Linear interpolant for non-uniform data. boost::barycentric_interpolation
    // yields nonsense and B-Splines require uniform data.

    auto P = [&](double k) {
        auto k1_i = discrete_P.lower_bound(k);
        auto k0_i = k1_i--;
        const double k0 = k0_i->first;
        const double P0 = k0_i->second;
        const double k1 = k1_i->first;
        const double P1 = k1_i->second;

        return P0 + (k - k0) * (P1 - P0) / (k1 - k0);
    };

    // Store fourier representation of delta in state.psis (complex vector)
    auto delta_k = subvector(state.psi, 0, N / 2 + 1);
    auto k = subvector(state.V, 0, N / 2 + 1);
    std::iota(k.begin(), k.end(), 0);

    // physical box size in Mpc h-1
    const double L_phys = cosmo.x_of_chi(L) / 1e6;
    k *= 2 * M_PI / L_phys;

    // Independent number generators for modulus and phase.
    // Seed chosen at random.
    auto uniform_rayleigh =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(seed));
    auto uniform_phase =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(seed));

    const double a_init = cosmo.a_of_tau(0);
    const double G = cosmo.Dplus(1) / cosmo.Dplus(a_init);

    auto sigma = [&](double k) {
        return sqrt(k * k / (4 * M_PI * L_phys * G * G) * P(k));
    };

    delta_k = map(delta_k, k, [&](std::complex<double> d, double k) {
        // Exact inversion yields transformation formulas for modulus and phase.
        return std::polar(sigma(k) * sqrt(-2 * log(uniform_rayleigh())),
                          2 * M_PI * uniform_phase());
    });

    auto in = reinterpret_cast<fftw_complex*>(delta_k.data());
    // Store delta in V for convenience (real vector vs complex vector);
    auto c2r = fftw_plan_dft_c2r_1d(N, in, state.V.data(), FFTW_ESTIMATE);
    fftw_execute(c2r);
    fftw_destroy_plan(c2r);
    */
}
