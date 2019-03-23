#include "ic.h"
#include <blaze/Blaze.h>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include "debugging.h"

ICGenerator::ICGenerator(const Parameters& param)
    : type(param.ic),
      N(param.N),
      dx(param.dx),
      L(param.L),
      M(param.M),
      rel_threshold(param.ev_thr),
      source_name(param.ic_source_file),
      ic_file(source_name) {
    // Altpugh the input file format is generic for rho and powerspectrum
    // type initial conditions, we leave the exact way of parsing the data to
    // the generator routines to allow for specialized malloc's such as done by
    // FFTW3. Here, we only discard the file header and count lines
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
    if (type == Parameters::IC::Density)
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

    // Compute FFT to obtain coefficients for the expansion of rho in
    // harmonic base functions
    auto plan = fftw_plan_dft_r2c_1d(
        data_N, rho.data(), reinterpret_cast<fftw_complex*>(rho_fft.data()),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // We only care about dominant modes in the statistic mixture, hence we
    // sort all eigenvalues, the modulus of the fourier coefficients, by
    // magnitude. In fact, we also need information about the mode number n to
    // sample from the correct harmonics later on. Therefore, we compute the
    // permutation index array that yields a sorted coefficient list.
    DynamicVector<int> mode(rho_fft.size() - 1);
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
    // No. of wavefunction different from M_true as we still have to take the DC
    // signal into account modeling the homogeneous background
    M_true += 1;

    // No .of wavefuntions with the same eigenvalue sign
    double M_half = (M_true - 1) / 2.0;

    std::cout << INFOTAG("Construct " << M_true << " wavefunctions")
              << std::endl;

    // At this point the matrix sizes for psi and V are clear. Set them.
    state.M = M_true;
    state.psis.resize(M_true, N);
    state.Vs.resize(M_true, N);
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
    state.lambda[M_true - 1] = 1.0 / data_N * rho_fft[0].real();

    // DC mode wavefunction is trivial
    row(state.psis, M_true - 1) = 1.0;

    // Equidistant x grid
    DynamicVector<double, rowVector> x(N);
    for (int i = 0; i < N; ++i) x[i] = dx * i;

    // Cosinus coefficients
    DiagonalMatrix<DynamicMatrix<double>> alpha(M_half, M_half);
    // Sinus coefficients
    DiagonalMatrix<DynamicMatrix<double>> beta(M_half, M_half);

    for (int i = 0; i < 2; ++i) {
        auto lambda = subvector(state.lambda, M_half * i, M_half);
        auto psi = submatrix(state.psis, M_half * i, 0, M_half, N);
        diagonal(alpha) = 2.0 / data_N * blaze::real(rho_fft_trunc) + lambda;
        diagonal(beta) = -2.0 / data_N * blaze::imag(rho_fft_trunc);

        // Normalization - decldiag() ???
        auto N = blaze::invsqrt(alpha * alpha + beta * beta);

        // Construct initial wave functions
        psi = N * (alpha * blaze::cos(M_PI / L * mode_trunc * x) +
                   beta * blaze::sin(M_PI / L * mode_trunc * x));
    }

    for (int m = 0; m < M_true - 1; m += 2) {
        // both solution have the same wave number
        int n = mode[m / 2];
        // Extract cos() and sin() coefficient from exp() coefficient
        double a = 2.0 / data_N * rho_fft[n].real();
        double b = -2.0 / data_N * rho_fft[n].imag();

        double lambda = eigenvalues[n];
        double normal = 1.0 / std::sqrt((a + lambda) * (a + lambda) + b * b);
        // Append (+)-solution
        for (int i = 0; i < N; ++i) {
            state.psis.push_back(
                normal * ((a + lambda) * std::cos(M_PI * n * dx * i / L) +
                          b * std::sin(M_PI * n * dx * i / L)));
        }
        state.lambda.push_back(lambda);

        lambda *= -1.0;
        normal = 1.0 / std::sqrt((a + lambda) * (a + lambda) + b * b);
        // Append (-)-solution
        for (int i = 0; i < N; ++i) {
            state.psis.push_back(
                normal * ((a + lambda) * std::cos(M_PI * n * dx * i / L) +
                          b * std::sin(M_PI * n * dx * i / L)));
        }
        state.lambda.push_back(lambda);
    }
    // Add constant wavefunction for background
    for (int i = 0; i < N; ++i) {
        state.psis.push_back(1.0);
    }
    state.lambda.push_back(1.0 / data_N * rho_fft[0].real());

    // Test
    for (int i = 0; i < M_true; ++i)
        for (int n = 0; n < N; ++n)
            state.Vs[n] += state.lambda[i] * norm(state.psis[n + i * N]);
}
void ICGenerator::psi_from_power(SimState& state) const {}
