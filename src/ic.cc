#include "ic.h"
#include "blaze/math/Band.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/Elements.h"
#include "blaze/math/Row.h"
#include "blaze/math/Submatrix.h"
#include "interfaces.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

// Stores a two column file into vectors col1 & col2. Copy process stops if end
// of vector is reached. Caller must ensure that stream is valid for at least
// 2*N reads from it.
// TODO: Regex support to generalize readout
template <typename T, bool TF>
void fill_from_file(std::istream& s, blaze::DynamicVector<T, TF>& col1,
                    blaze::DynamicVector<T, TF>& col2) {
    assert(col1.size() == col2.size());
    const size_t N = col1.size();
    for (size_t i = 0; i < N; ++i) {
        T val1, val2;
        s >> col1[i] >> col2[i];
    }
}

ICGenerator::ICGenerator(const Parameters& p)
    : type{static_cast<ICType>(p["Initial Conditions"]["ic_type"].get<int>())},
      N{p["Simulation"]["N"].get<int>()},
      data_N{0},
      L{p["Simulation"]["L"].get<double>()},
      dx{L / N},
      M{p["Initial Conditions"]["M"].get<int>()},
      rel_threshold{p["Initial Conditions"]["ev_threshold"].get<double>()},
      ic_file{p["Initial Conditions"]["source_file"].get<std::string>()},
      pot_file{p["Initial Conditions"]["potential_file"].get<std::string>()},
      potential{nullptr} {
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');
    data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                        std::istreambuf_iterator<char>(), '\n');

    int pot_N;
    if (pot_file) {
        pot_file.seekg(0);
        pot_file.ignore(std::numeric_limits<int>::max(), '\n');

        pot_N = std::count(std::istreambuf_iterator<char>(pot_file),
                           std::istreambuf_iterator<char>(), '\n');

        pot_file.seekg(0);
        pot_file.ignore(std::numeric_limits<int>::max(), '\n');
    } else {
        potential = PotentialMethod::make(
            p["Simulation"]["potential"].get<std::string>(), p);
        pot_N = N;
    }

    switch (type) {
        case ICType::External:
            if (data_N % M != 0) {
                std::cout << ERRORTAG(
                                 "Number of lines in source file is not a "
                                 "multiple of M.")
                          << std::endl;
                exit(1);
            }
            break;
        case ICType::Density:
            if (data_N != N) {
                std::cout << ERRORTAG(
                                 "Number of lines in source file differs from "
                                 "parameter N.")
                          << std::endl;
                exit(1);
            }
            break;
        default:
            if (N != pot_N) {
                std::cout
                    << ERRORTAG(
                           "Number of lines in potential file differs from N.")
                    << std::endl;
                exit(1);
            }
    };

    ic_file.seekg(0);
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');
}

void ICGenerator::generate(SimState& state, Parameters& param) const {
    switch (type) {
        case ICType::External:
            psi_from_file(state);
            break;
        case ICType::Powerspectrum:
            psi_from_power(state);
            break;
        case ICType::Density:
            psi_from_rho(state);
    }

    param["Initial Conditions"]["M"] = state.M;

    if (pot_file) {
        fill_from_file(pot_file, state.V, state.V);
    } else {
        potential->solve(state);
    }
}

void ICGenerator::psi_from_file(SimState& state) const {
    using namespace blaze;
    state.M = M;
    state.psis.resize(N, M);
    state.V.resize(N);
    state.lambda.resize(M);

    DynamicVector<double> real(N);
    DynamicVector<double> imag(N);

    for (int m = 0; m < M; ++m) {
        state.lambda[m] = 1.0;
        auto psi = column(state.psis, m);
        fill_from_file(ic_file, real, imag);
        psi = map(real, imag,
                  [](double r, double i) { return std::complex(r, i); });
    }
}

// TODO implement power spectrum initial conditions.
void ICGenerator::psi_from_power(SimState& state) const {}

void ICGenerator::psi_from_rho(SimState& state) const {
    using namespace blaze;

    // Read in density data from file
    DynamicVector<double> rho(N);
    fill_from_file(ic_file, rho, rho);

    // Compute FFT to obtain coefficients for the expansion of rho in
    // harmonic base functions
    DynamicVector<std::complex<double>> rho_fft(N / 2 + 1);
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
    DynamicVector<int, rowVector> mode(rho_fft.size() - 1);
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

    // Constant wvaefunction for background
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
