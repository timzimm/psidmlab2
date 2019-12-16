#include "ic.h"
#include "cosmology.h"
#include "fftw.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>

using namespace blaze;

ICGenerator::ICGenerator(const Parameters& p)
    : type{static_cast<ICType>(p["Initial Conditions"]["ic_type"].get<int>())},
      N{p["Simulation"]["N"].get<int>()},
      data_N{0},
      L{p["Simulation"]["L"].get<double>()},
      dx{L / N},
      seed{p["Initial Conditions"]["seed"].get<int>()},
      compute_velocity{p["Initial Conditions"]["compute_velocity"].get<bool>()},
      ic_file{p["Initial Conditions"]["source_file"].get<std::string>()},
      pot{nullptr} {
    data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                        std::istreambuf_iterator<char>(), '\n');
    if (compute_velocity) {
        pot = Interaction::make("Poisson::FFT", p);
    }

    switch (type) {
        case ICType::ExternalRealImag... ICType::ExternalModulusPhase:
            if (data_N != N) {
                std::cout << ERRORTAG("#lines in source_file(" << data_N
                                                               << ") != N")
                          << std::endl;
                exit(1);
            }
        default:
            break;
    };

    ic_file.seekg(0);
}

void ICGenerator::generate(SimState& state, const Cosmology& cosmo) const {
    // modulus initialization
    switch (type) {
        case ICType::ExternalRealImag:
            std::cout << INFOTAG("Load Re(psi0) and Im(psi0) from file")
                      << std::flush;
            real_imag_from_file(state);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::ExternalModulusPhase:
            std::cout << INFOTAG("Load psi0 from file") << std::flush;
            modulus_phase_from_file(state);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::Powerspectrum:
            std::cout << INFOTAG("Generate delta0 from P(k)") << std::flush;
            delta_from_power(state, cosmo);
            std::cout << " ... done" << std::endl;
            break;
    }

    // Velocity Initialization
    if (type == ICType::Powerspectrum) {
        // state.V holds delta at this point. But we will use it in an in-place
        // manner. Hence we define...
        const auto& delta = state.V;
        auto& phase = state.V;
        auto& psi = state.psi;

        // Set psis modulus before we override delta. For cold conditions that's
        // all that is left to do.
        psi = sqrt(1.0 + delta);

        // Cosmological initial velocity field
        if (compute_velocity) {
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

void ICGenerator::real_imag_from_file(SimState& state) const {
    state.psi.resize(N);
    state.V.resize(N);

    DynamicVector<double> mod(N);
    DynamicVector<double> phase(N);

    fill_from_file(ic_file, mod, phase);
    state.psi =
        map(mod, phase, [](double m, double p) { return std::polar(m, p); });
}

void ICGenerator::modulus_phase_from_file(SimState& state) const {
    state.psi.resize(N);
    state.V.resize(N);

    DynamicVector<double> re(N);
    DynamicVector<double> im(N);

    fill_from_file(ic_file, re, im);
    state.psi =
        map(re, im, [](double r, double i) { return std::complex(r, i); });
}

void ICGenerator::delta_from_power(SimState& state,
                                   const Cosmology& cosmo) const {
    state.psi.resize(N);
    state.V.resize(N);

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
    fftw_plan_ptr c2r(
        fftw_plan_dft_c2r_1d(N, in, state.V.data(), FFTW_ESTIMATE));
    fftw_execute(c2r.get());
    fftw_destroy_plan(c2r.get());
}

ICGenerator::~ICGenerator() = default;
