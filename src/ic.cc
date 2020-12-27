#include "ic.h"
#include "config.h"
#include "cosmology.h"
#include "fftw.h"
#include "hdf5_file.h"
#include "interaction/periodic_convolution.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "parameters_fwd.h"
#include "state.h"

#include <algorithm>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <random>

using namespace blaze;

ICGenerator::ICGenerator(const Parameters &p)
    : type{static_cast<ICType>(p["Initial Conditions"]["ic_type"])},
      data_N(0), box{p}, seed(0),
      compute_velocity(false), filename{
                                   p["Initial Conditions"]["source_file"]} {
    if (type != ICType::PreviousSimulation) {
        ic_file = std::ifstream(filename);
        data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                            std::istreambuf_iterator<char>(), '\n');
        ic_file.seekg(0);
    }
    if (type < ICType::Powerspectrum) {
        if (data_N != box.N) {
            std::cerr << ERRORTAG("#lines in source_file(" << data_N
                                                           << ") != N")
                      << std::endl;
            exit(1);
        }
    }
    if (type == ICType::Powerspectrum) {
        if (p["Domain"]["physical_units"] == false) {
            std::cerr << ERRORTAG(
                             "psi0 from power spectrum requires physical units")
                      << std::endl;
            exit(1);
        }
        compute_velocity = p["Initial Conditions"]["compute_velocity"];
        seed = p["Initial Conditions"]["seed"];

    } else if (type == ICType::PreviousSimulation) {
        auto isHDF5 = H5Fis_hdf5(filename.c_str());
        if (isHDF5 <= 0) {
            std::cerr << ERRORTAG("No valid HDF5 file found") << std::endl;
            exit(1);
        }
        HDF5File file(filename, HDF5File::Access::Read);

        std::vector<std::string> paths = file.ls("/WaveFunction/");
        if (paths.size() == 0) {
            std::cerr << ERRORTAG("No wavefunction found.") << std::endl;
            exit(1);
        }
        for (auto &path : paths) {
            psis.push_back(path);
            auto tau = file.read_scalar_attribute<double>(path, "tau");
            known_tau.push_back(tau);
        }
    }
}

void ICGenerator::generate(SimState &state, const Cosmology &cosmo,
                           const double tau) const {
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
    case ICType::PreviousSimulation:
        std::cout << INFOTAG("Load psi from state") << std::flush;
        psi_from_state(state, cosmo, tau);
    }

    // Velocity Initialization
    if (type == ICType::Powerspectrum) {
        // state.V holds delta at this point. But we will use it in an
        // in-place manner. Hence we define...
        auto &delta = state.V;
        auto &phase = state.V;
        auto &psi = state.psi;

        // Set psis modulus before we override delta. For cold conditions
        // that's all that is left to do.
        psi = blaze::sqrt(1.0 + delta);

        // Cosmological initial velocity field
        if (compute_velocity) {
            const Parameters p(
                {{"Domain",
                  {{"N", box.N}, {"L", box.L}, {"physical_units", false}}},
                 {"Simulation", {{"interaction", {{"type", 0}}}}}});
            const double a_init = cosmo.a_of_tau(0);
            const double prefactor =
                -std::sqrt(2.0 / 3 * a_init / cosmo.omega_m(a_init));

            std::cout << INFOTAG("Initial Velocity Field from Poisson @ a = ")
                      << a_init << std::endl;

            auto pot = Interaction::make("PeriodicConvolution", p, state);
            pot->solve(phase, prefactor * delta);

            // Madelung Representation.
            psi *= exp(std::complex<double>(0, 1) * phase);
        }
    }
}

// Populate SimState with stored wavefunction if the requested tau is greater
// than the stored tau.
void ICGenerator::psi_from_state(SimState &state, const Cosmology &cosmo,
                                 const double tau) const {
    static bool quick_exit = false;
    if (quick_exit) {
        return;
    }

    HDF5File file(filename, HDF5File::Access::Read);
    // Find next stored tau greater or equal to the requested tau - epsilon
    auto tau_up_iter =
        std::lower_bound(known_tau.begin(), known_tau.end(),
                         tau - PSIDMLAB_MINIMUM_INTEGRATION_STEP);

    double next_tau;
    std::string path_to_next_psi;
    // (i) Simulation Resumed:
    // Requested tau is beyond the last stored psi. In this case we
    // load the last state and make sure we exit immediatly if this function
    // will be called again.
    if (tau_up_iter == known_tau.end()) {
        next_tau = known_tau.back();
        path_to_next_psi = psis.back();
        quick_exit = true;
        // (ii) Simulation Refinement:
        // Requested tau lies within the time grid of known states. In this case
        // we load the state with largest state.tau <= tau to resume the
        // integration from.
    } else {
        double tau_up = *tau_up_iter;
        size_t idx_to_tau_up = std::distance(known_tau.begin(), tau_up_iter);
        // Load upper bound if difference is so small (but possibly negative)
        // that we would not integrate anyways. This happens if new
        // observables are computed at values of tau where the wavefunction is
        // already known
        if (std::abs(tau_up - tau) < PSIDMLAB_MINIMUM_INTEGRATION_STEP) {
            next_tau = tau_up;
            path_to_next_psi = psis[idx_to_tau_up];
            // Load lower bound if observables need to be stored at
            // values of tau at which no psi was previously stored.
        } else {
            double tau_low = *(tau_up_iter - 1);
            next_tau = tau_low;
            path_to_next_psi = psis[idx_to_tau_up - 1];
        }
    }
    state.tau = next_tau;
    state.n = file.read_scalar_attribute<unsigned int>(path_to_next_psi, "n");

    auto psi_matrix = file.read_matrix(path_to_next_psi);
    // Populate simualation state.
    state.psi =
        blaze::map(blaze::column(psi_matrix, 0), blaze::column(psi_matrix, 1),
                   [](const double r, const double i) {
                       return std::complex<double>(r, i);
                   });
    state.V.resize(box.N);
    std::cout << " @ z = " << Cosmology::z_of_a(cosmo.a_of_tau(state.tau))
              << ", t = " << state.tau << " ... done" << std::endl;
}

void ICGenerator::real_imag_from_file(SimState &state) const {
    fill_from_file(ic_file, state.psi);
}

void ICGenerator::modulus_phase_from_file(SimState &state) const {
    fill_from_file(ic_file, state.psi);
    state.psi = map(state.psi, [](std::complex<double> p) {
        return std::polar(p.real(), p.imag());
    });
    if (isnan(state.psi)) {
        std::cerr << ERRORTAG("Not a valid polar representation") << std::endl;
        exit(1);
    }
}

void ICGenerator::delta_from_power(SimState &state,
                                   const Cosmology &cosmo) const {
    using namespace boost::math::interpolators;
    // First line is a header. Ignore it.
    ic_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    size_t start = ic_file.tellg();
    int Nk = data_N - 1;

    std::vector<double> Pk(Nk);
    fill_from_file(ic_file, Pk, Pk);

    ic_file.seekg(start);
    double k0, k1;
    ic_file >> k0;
    ic_file >> k1 >> k1;
    const double logk0 = std::log10(k0);
    const double delta_logk = std::log10(k1) - logk0;

    // Cubic interpolation of 3D CDM spectrum in logk
    // CAMB returns P(k/h) in units of h^{-3} * Mpc^3
    // P_3d returns P(k) in units of Mpc^3
    // Note that both arguments have unit Mpc^{-1}
    auto P_3d = [P = cardinal_cubic_b_spline<double>(Pk.begin(), Pk.end(),
                                                     logk0, delta_logk),
                 kmin_div_h = k0,
                 kmax_div_h = std::pow(10, logk0 + (Nk - 1) * delta_logk),
                 h = box.hubble](double k) {
        const double k_div_h = k / h;
        const double h3 = 1.0 / (h * h * h);
        // Power law extrapolation before kmin. Assume P(k) = A*k
        if (k_div_h < kmin_div_h) {
            return h3 * P(std::log10(kmin_div_h)) / (kmin_div_h * h) * k;
        }
        // Power law extrapolation past kmax. Assume P(k) = A*k^-3
        if (k_div_h > kmax_div_h) {
            return h3 * P(std::log10(kmax_div_h)) *
                   std::pow(kmax_div_h * h, 3) * std::pow(k, -3);
        }
        return h3 * P(std::log10(k_div_h));
    };

    // Store fourier representation of delta in state.psis (complex vector)
    auto delta_k = subvector(state.psi, 0, box.N / 2 + 1);
    double L_Mpc = box.L_phys;
    // Modes in units of Mpc^-1
    auto k =
        blaze::linspace(box.N / 2 + 1, 0.0, 2 * M_PI / L_Mpc * (box.N / 2));

    // Independent number generators for modulus and phase.
    // Seed chosen at random.
    auto uniform1 =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(seed));
    auto uniform2 =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(seed));
    const double a_init = cosmo.a_of_tau(0);
    const double G = cosmo.D(a_init) / cosmo.D(1);

    auto sigma = [&](double k) {
        double kJeq = 9 * std::sqrt(box.m22);
        double x = 1.61 * std::pow(box.m22, 1.0 / 18) * k / kJeq;
        double TFDM = std::cos(x * x * x) / (1 + std::pow(x, 8));
        return sqrt(k * k * G * G * TFDM * TFDM * L_Mpc / (4 * M_PI) * P_3d(k));
    };

    delta_k = map(k, [&](double k) {
        // Exact inversion yields transformation formulas for modulus and
        // phase.
        return std::polar(sigma(k) * std::sqrt(-2 * std::log(uniform1())),
                          2 * M_PI * uniform2());
    });
    if (box.N % 2 == 0) {
        delta_k[box.N / 2] = std::sqrt(2) * sigma(k[box.N / 2]) *
                             std::sqrt(-2 * std::log(uniform1())) *
                             std::cos(2 * M_PI * uniform2());
    }

    // Store delta in state.V for convenience (real vector)
    auto &delta = state.V;
    auto in = reinterpret_cast<fftw_complex *>(delta_k.data());

    fftw_plan_ptr c2r(
        fftw_plan_dft_c2r_1d(box.N, in, delta.data(), FFTW_ESTIMATE));
    fftw_execute(c2r.get());
    // Note that delta_k has dimensions of L. It is related to the DFT
    // coeffcient via delta_k = delta_k_DFT * L_Mpc
    delta /= L_Mpc;
}
