#include "observables.h"
#include <complex>
#include <iostream>
#include "blaze/math/UniformVector.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include "blaze/math/Elements.h"
#include "blaze/math/Rows.h"
#include "blaze/math/Submatrix.h"

namespace Observable {

// TODO The initializer list is too messy
DensityContrast::DensityContrast(const Parameters& p)
    : sigma_x(p["Analysis"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      ws(husimi ? p["Analysis"]["linear_convolution"].get<bool>() : 0,
         husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(N_kernel) {
    // Kernel construction in x-space
    for (int i = 0; i < N_kernel; ++i) {
        gaussian_kernel[i] = (i - N_kernel / 2) * dx;
    }
    gaussian_kernel =
        exp(-0.5 / (sigma_x * sigma_x) * gaussian_kernel * gaussian_kernel);
    gaussian_kernel /= blaze::sum(gaussian_kernel);
}

blaze::DynamicMatrix<double, blaze::columnMajor> DensityContrast::compute(
    const SimState& state) {
    auto psi2 = blaze::real(state.psis % state.psis);
    blaze::UniformVector<double> one(psi2.rows());
    blaze::DynamicVector<double> delta = psi2 * state.lambda - one;

    if (husimi) {
        discrete_convolution(ws, gaussian_kernel, delta);
        delta = blaze::subvector(ws.signal_padded, ws.N_offset, N);
    }

    return blaze::expand(delta, 1);
}

PhaseSpaceDistribution::PhaseSpaceDistribution(const Parameters& p)
    : sigma_x(p["Analysis"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      ws(husimi ? p["Analysis"]["linear_convolution"].get<bool>() : 0,
         husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(N_kernel) {}

blaze::DynamicMatrix<double, blaze::columnMajor>
PhaseSpaceDistribution::compute(const SimState& state) {
    blaze::DynamicVector<std::complex<double>> chirp(N);
    using namespace blaze;

    // Load chirp data
    std::ifstream input("/tmp/chirp/chirp.txt");
    std::copy(std::istream_iterator<std::complex<double>>(input),
              std::istream_iterator<std::complex<double>>(), chirp.begin());

    auto& signal = chirp;

    DynamicMatrix<double, columnMajor> W(2 * N, 2 * N);

    // We only compute the first quadrant of W. All other quadrants follow by
    // symmetry properties of W
    auto W_00 = blaze::submatrix(W, 0, 0, N, N);

    DynamicVector<std::complex<double>, blaze::rowVector> strip(N);

    auto in_out = reinterpret_cast<fftw_complex*>(strip.data());

    fftw_plan one =
        fftw_plan_dft_1d(N, in_out, in_out, FFTW_FORWARD, FFTW_ESTIMATE);

    // First we create the "lag matrix" given by
    //
    // W(k, 2*n) = 1/(2N) * signal[(k+n)modN]signal[(n-k)modN]^*
    // W(k, 2*n+1) = 1/(2N) * signal[(k+n)modN]signal[(n-k+1)modN]^*

    for (int n = 0; 2 * n < N; ++n) {
        auto sig_npk =
            elements(signal, [&](size_t k) { return (n + k) % N; }, N);
        auto sig_nmk =
            elements(signal, [&](size_t k) { return (n - k) % N; }, N);
        auto sig_nmk1 =
            elements(signal, [&](size_t k) { return (n - k + 1) % N; }, N);

        strip = 1.0 / (2 * N) * sig_npk * conj(sig_nmk);
        fftw_execute(one);
        column(W_00, 2 * n) = real(strip);
        strip = 1.0 / (2 * N) * sig_npk * conj(sig_nmk1);
        fftw_execute(one);
        for (int m = 0; m < N; ++m)
            strip[m] *= std::exp(std::complex<double>(0, 1) * (M_PI * m / N));
        column(W_00, 2 * n + 1) = real(strip);
    }

    // W(k,n+N) = (-1)^k W(k,n)
    auto W_01 = submatrix(W, 0, N, N, N) = W_00;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; 2 * k < N; ++k) {
            W_01(2 * k + 1, n) *= -1;
        }
    }

    // W(k+N,n) = (-1)^n W(k,n)
    auto W_10 = submatrix(W, N, 0, N, N) = W_00;
    for (int n = 0; 2 * n < N; ++n) {
        for (int k = 0; 2 * k < N; ++k) {
            W_10(k, 2 * n + 1) *= -1;
        }
    }

    // W(k+N,n+N) = (-1)^(k+n+N) W(k,n)
    auto W_11 = submatrix(W, N, N, N, N) = W_00;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            W_11(k, n) *= std::pow(-1, k + n + N);
        }
    }

    return W;
}

}  // namespace Observable

