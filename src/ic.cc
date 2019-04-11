#include "ic.h"
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include "blaze/math/Band.h"
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/Elements.h"
#include "blaze/math/Row.h"
#include "blaze/math/Submatrix.h"
#include "common.h"
#include "interfaces.h"
#include "logging.h"

ICGenerator::ICGenerator(const Parameters& param)
    : type(param.ic),
      N(param.N),
      dx(param.dx),
      L(param.L),
      M(param.M),
      rel_threshold(param.ev_thr),
      source_name(param.ic_source_file),
      ic_file(source_name),
      potential(PotentialMethod::make("Poisson::FD", param)) {
    // Altough the input file format is generic for rho and powerspectrum
    // type initial conditions, we leave the exact way of parsing the data
    // to the generator routines to allow for specialized malloc's.
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');

    // Number of lines
    data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                        std::istreambuf_iterator<char>(), '\n');
    data_N = (data_N + 1);
    ic_file.seekg(0);
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');
}

// Dispatches the correct generator routine depending on
// Parameters::IC parameter at runtime
void ICGenerator::generate(SimState& state) const {
    if (type == ICType::Density)
        psi_from_rho(state);
    else
        psi_from_power(state);
}

void ICGenerator::psi_from_rho(SimState& state) const {
    using namespace blaze;

    // Read in density data from file
    DynamicVector<double> rho(data_N);
    DynamicVector<std::complex<double>> rho_fft(data_N / 2 + 1);
    for (int i = 0; i < data_N; ++i) {
        ic_file >> rho[i] >> rho[i];
        if (!ic_file) break;
    }

    // TODO Move into constructor
    // Compute FFT to obtain coefficients for the expansion of rho in
    // harmonic base functions
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
        std::cout << INFOTAG("Determine No. of wavefunctions by fiven M")
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

    for (int i = 0; i < 2; ++i) {
        auto lambda = subvector(state.lambda, M_half * i, M_half);
        auto psi = submatrix(state.psis, 0, M_half * i, N, M_half);
        diagonal(alpha) = 2.0 / data_N * real(rho_fft_trunc) + lambda;
        diagonal(beta) = -2.0 / data_N * imag(rho_fft_trunc);

        // Normalization
        auto N = decldiag(invsqrt(alpha * alpha + beta * beta));

        // Construct initial wave functions
        psi = (cos(M_PI / L * x * mode_trunc) * alpha +
               sin(M_PI / L * x * mode_trunc) * beta) *
              N;
    }

    // Compute potential
    potential->solve(state);
}

// TODO implement power spectrum initial conditions.
void ICGenerator::psi_from_power(SimState& state) const {}

ICGenerator::~ICGenerator() = default;
