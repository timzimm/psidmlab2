#include "ic.h"
#include "cosmology.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

#include <blaze/math/Band.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/Elements.h>
#include <blaze/math/Row.h>
#include <blaze/math/Submatrix.h>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <numeric>

// Reads columns of stream s into the supplied vectors.
template <typename... Ts>
void fill_from_file(std::istream& s, Ts&... vecs) {
    const int N = std::count(std::istreambuf_iterator<char>(s),
                             std::istreambuf_iterator<char>(), '\n');
    s.seekg(0);
    for (int i = 0; i < N; ++i) {
        (s >> ... >> vecs[i]);
    }
}

ICGenerator::ICGenerator(const Parameters& p)
    : type{static_cast<ICType>(p["Initial Conditions"]["ic_type"].get<int>())},
      N{p["Simulation"]["N"].get<int>()},
      data_N{0},
      L{p["Simulation"]["L"].get<double>()},
      dx{L / N},
      M{(type == ICType::Experimental) ? p["Initial Conditions"]["M"].get<int>()
                                       : 1},
      rel_threshold{p["Initial Conditions"]["ev_threshold"].get<double>()},
      compute_velocity{p["Initial Conditions"]["compute_velocity"].get<bool>()},
      ic_file{p["Initial Conditions"]["source_file"].get<std::string>()},
      pot_file{p["Initial Conditions"]["potential_file"].get<std::string>()},
      poisson{nullptr} {
    data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                        std::istreambuf_iterator<char>(), '\n');
    int pot_N;
    if (pot_file) {
        pot_file.seekg(0);

        pot_N = std::count(std::istreambuf_iterator<char>(pot_file),
                           std::istreambuf_iterator<char>(), '\n');

        pot_file.seekg(0);
    } else {
        poisson = PotentialMethod::make(
            p["Simulation"]["potential"].get<std::string>(), p);
        pot_N = N;
    }
    if (compute_velocity && !poisson) {
        poisson = PotentialMethod::make("Poisson::FFT", p);
    }

    switch (type) {
        case ICType::External:
            if (data_N != N * M) {
                std::cout << ERRORTAG("#lines in source_file(" << data_N
                                                               << ") != N*M")
                          << std::endl;
                exit(1);
            }
        case ICType::Experimental:
            if (data_N != N) {
                std::cout << ERRORTAG("#lines in source_file(" << data_N
                                                               << ") != N")
                          << std::endl;
                exit(1);
            }
        default:
            if (N != pot_N) {
                std::cout << ERRORTAG("#lines in pot_file(" << pot_N
                                                            << ") != N")
                          << std::endl;
                exit(1);
            }
    };

    ic_file.seekg(0);
}

void ICGenerator::generate(SimState& state, const Cosmology& cosmo) const {
    using namespace blaze;
    // delta initialization
    switch (type) {
        case ICType::External:
            std::cout << INFOTAG("Load delta0 from file") << std::endl;
            delta_from_file(state);
            break;
        case ICType::Powerspectrum:
            std::cout << INFOTAG("Generate delta0 from P(k)") << std::endl;
            delta_from_power(state);
            break;
        case ICType::Experimental:
            std::cout << INFOTAG("Generate Psi0 from rho(x)") << std::endl;
            delta_from_evp(state);
    }

    // Velocity Initialization - currently for pure states only
    if (type != ICType::Experimental) {
        // state.V holds delta at this point. But we will use it in an in-place
        // manner. Hence we define...
        const auto& delta = state.V;
        auto& phase = state.V;
        auto psi = column(state.psis, 0);

        // Set psis modulus before we override delta. For cold conditions that's
        // all that is left to do.
        psi = sqrt(1.0 + delta);

        if (compute_velocity) {
            const double a_init = cosmo.a_of_tau(0);
            const double prefactor =
                -std::sqrt(2.0 / 3 * a_init / cosmo.omega_m(a_init));

            std::cout << INFOTAG("Initial Velocity Field from Poisson @ a = ")
                      << a_init << std::endl;

            poisson->solve(phase, prefactor * delta);

            // Madelung Representation.
            psi *= exp(std::complex<double>(0, 1) * phase);
        }
    }

    // Potential Initialization
    if (pot_file) {
        std::cout << INFOTAG("Load potential from file.") << std::endl;
        fill_from_file(pot_file, state.V);
    } else {
        std::cout << INFOTAG("Generate Potential from Poisson") << std::endl;
        poisson->solve(state);
    }
}

// currently for pure states only (M=1)
void ICGenerator::delta_from_file(SimState& state) const {
    state.M = M;
    state.psis.resize(N, M);
    state.V.resize(N);
    state.lambda.resize(M);

    state.lambda[0] = 1.0;
    // Store delta in V for convenience (real vector vs complex vector);
    fill_from_file(ic_file, state.V);
}

// TODO implement power spectrum initial conditions.
// currently for pure states only (M=1)
void ICGenerator::delta_from_power(SimState& state) const {}

void ICGenerator::delta_from_evp(SimState& state) const {
    using namespace blaze;

    // Read in density data from file
    DynamicVector<double> rho(data_N);
    fill_from_file(ic_file, rho);

    // Compute FFT to obtain coefficients for the expansion of rho in
    // harmonic base functions
    DynamicVector<std::complex<double>> rho_fft(data_N / 2 + 1);
    auto plan = fftw_plan_dft_r2c_1d(
        data_N, rho.data(), reinterpret_cast<fftw_complex*>(rho_fft.data()),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // We only care about dominant modes in the statistic mixture, hence we
    // sort all eigenvalues, the modulus of the fourier coefficients, by
    // magnitude. In fact, we also need information about the mode number n
    // to sample from the correct harmonics later on. Therefore, we compute
    // the permutation index array that yields a sorted coefficient list.
    DynamicVector<int, rowVector> mode(data_N - 1);
    std::iota(mode.begin(), mode.end(), 1);

    // Construct the permuatation/mode-number vector
    std::sort(mode.begin(), mode.end(), [&](const int& a, const int& b) {
        return norm(rho_fft[a]) > norm(rho_fft[b]);
    });

    // Sort the actual coefficient vector according to the permutation
    auto rho_fft_sorted = elements(rho_fft, mode.data(), mode.size());

    // Determine dominant modes by threshold or max number of wavefunctions
    int M_true;
    // (i) rel_threshold <= 0 -> min(M, 2 * #fourier_modes)
    if (rel_threshold <= 0) {
        std::cout << INFOTAG("Determine No. of wavefunctions by given M")
                  << std::endl;
        M_true = std::min(M, 2 * static_cast<int>(rho_fft_sorted.size()));
        M_true += (M_true % 2 == 0) ? 0 : 1;
    }

    // (ii) rel_threshold > 0 -> min(M s.t. all ev >
    // ev_threshold,#fourier_modes)
    else {
        double abs_thr = rel_threshold * std::norm(rho_fft_sorted[0]);
        std::cout << INFOTAG(
                         "Determine No. of wavefunctions by given threshold")
                  << std::endl;
        std::cout << INFOTAG("Absolut coefficient threshold: " << abs_thr)
                  << std::endl;

        // min(..) is already taken care of.
        M_true = std::count_if(rho_fft_sorted.begin(), rho_fft_sorted.end(),
                               [&](const std::complex<double>& c) {
                                   return abs_thr < std::norm(c);
                               });

        M_true *= 2;
    }

    // No .of wavefuntions with the same eigenvalue sign
    int M_half = M_true / 2;

    // Constant wavefunction for background
    M_true += 1;

    std::cout << INFOTAG("Construct " << M_true << " wavefunctions")
              << std::endl;

    // At this point the matrix sizes for psi and V are clear. Set them.
    state.M = M_true;
    state.psis.resize(N, M_true);
    state.V.resize(N);
    state.lambda.resize(M_true);

    // Discard all irrelvant Fourier coefficients and modes
    auto rho_fft_trunc = subvector(rho_fft_sorted, 0, M_half);
    auto mode_trunc = subvector(mode, 0, M_half);

    // The rescaled modulus of the Fourier coefficients are the eigenvalues
    // of the integral operator eigenvalue problem to solve.
    auto lambda_plus = subvector(state.lambda, 0, M_half);
    auto lambda_minus = subvector(state.lambda, M_half, M_half);
    lambda_plus = 2.0 / data_N * abs(rho_fft_trunc);
    lambda_minus = -2.0 / data_N * abs(rho_fft_trunc);

    // Equidistant x grid
    DynamicVector<double, columnVector> x(N);
    for (int i = 0; i < N; ++i) x[i] = dx * i;

    // TODO switch to Sparse Matrix ?
    // Cosinus coefficients
    DiagonalMatrix<DynamicMatrix<double>> alpha(M_half, M_half);
    // Sinus coefficients
    DiagonalMatrix<DynamicMatrix<double>> beta(M_half, M_half);
    // Normalization
    DiagonalMatrix<DynamicMatrix<double>> normal(M_half, M_half);

    for (int i = 0; i < 2; ++i) {
        auto lambda = subvector(state.lambda, M_half * i, M_half);
        auto psi = submatrix(state.psis, 0, M_half * i, N, M_half);
        diagonal(alpha) = 2.0 / data_N * real(rho_fft_trunc) + lambda;
        diagonal(beta) = -2.0 / data_N * imag(rho_fft_trunc);

        // Normalization
        auto normsq = decldiag(alpha * alpha + beta * beta);
        diagonal(normal) = invsqrt(diagonal(normsq));

        // Construct initial wave functions
        psi = (cos(M_PI / L * x * mode_trunc) * alpha +
               sin(M_PI / L * x * mode_trunc) * beta) *
              normal;
    }

    // Constant Wavefunction
    column(state.psis, M_true - 1) = 1;
    state.lambda[M_true - 1] = 1;
}

ICGenerator::~ICGenerator() = default;
