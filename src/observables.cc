#include "observables.h"
#include <complex>
#include <iostream>
#include "blaze/math/Columns.h"
#include "blaze/math/Elements.h"
#include "blaze/math/Rows.h"
#include "blaze/math/Submatrix.h"
#include "blaze/math/UniformVector.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>

namespace Observable {

DensityContrast::DensityContrast(const Parameters& p)
    : sigma_x(p["Analysis"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Analysis"]["linear_convolution"].get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(0) {
    if (husimi) {
        gaussian_kernel.resize(N_kernel);
        // Kernel construction in x-space
        for (int i = 0; i < N_kernel; ++i) {
            gaussian_kernel[i] = (i - N_kernel / 2) * dx;
        }
        gaussian_kernel =
            exp(-0.5 / (sigma_x * sigma_x) * gaussian_kernel * gaussian_kernel);
        gaussian_kernel /= blaze::sum(gaussian_kernel);
    }
}

blaze::DynamicMatrix<double, blaze::columnMajor> DensityContrast::compute(
    const SimState& state) {
    auto psi2 = blaze::real(state.psis % state.psis);
    blaze::UniformVector<double> one(psi2.rows());
    blaze::DynamicVector<double> delta = psi2 * state.lambda - one;

    if (husimi) {
        discrete_convolution(ws, gaussian_kernel, delta);
        delta = blaze::subvector(ws.signal_padded, ws.N_offset, N);
        if (!linear) {
            // TODO Benchmark against element selection
            std::rotate(std::begin(delta), std::begin(delta) + N_kernel / 2,
                        std::end(delta));
        }
    }

    return blaze::expand(delta, 1);
}

PhaseSpaceDistribution::PhaseSpaceDistribution(const Parameters& p)
    : sigma_x(p["Analysis"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Analysis"]["linear_convolution"].get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      gaussian_kernel(N_kernel) {}

blaze::DynamicMatrix<double, blaze::columnMajor>
PhaseSpaceDistribution::compute(const SimState& state) {
    using namespace blaze;

    // Plan init is negligible compared to workload - do it here
    DynamicVector<std::complex<double>> lag(N);
    DynamicMatrix<double, columnMajor> W(2 * N, 2 * N);
    auto lag_ptr = reinterpret_cast<fftw_complex*>(lag.data());

    auto even_plan = fftw_plan_dft_c2r_1d(N, lag_ptr, W.data(), FFTW_ESTIMATE);
    auto odd_plan =
        fftw_plan_dft_1d(N, lag_ptr, lag_ptr, FFTW_FORWARD, FFTW_ESTIMATE);

    // We only compute the first quadrant of W. All other quadrants follow by
    // symmetry properties of W
    auto W_00 = submatrix(W, 0, 0, N, N);

    for (int m = 0; m < state.M; ++m) {
        auto psi = column(state.psis, m);
        // Create Wigner distribution column by column.

        // TODO: We could, in principle, first construct the entire lag matrix
        // and then transform it once via FFTWs advanced interface. This would
        // increase preformance (potentially significantly). The downside is
        // that we either have to store two matrices to do the out of place FFT,
        // or have one padded matrix for which we have to mess with C++ type
        // system to accomodate std::complex and doubles at once. I don't know
        // if this is really the right approach and worth the effort.

        for (int n = 0; 2 * n < N; ++n) {
            // First we create the "lag matrix" column given by
            // W(k, 2*n) = 1/(2N) * psi[(k+n)modN]psi[(n-k)modN]^*
            // W(k, 2*n+1) = 1/(2N) * psi[(k+n)modN]psi[(n-k+1)modN]^*
            // Only works for power of 2 N's !!!
            auto sig_npk =
                elements(psi, [&](size_t k) { return (n + k) % N; }, N);
            auto sig_nmk =
                elements(psi, [&](size_t k) { return (n - k) % N; }, N);
            auto sig_nmk1 =
                elements(psi, [&](size_t k) { return (n - k + 1) % N; }, N);

            // Even columns
            // conj switched due to sign convention of FFTW. Result is real so
            // no additional conj required.
            lag = state.lambda[m] / (2 * N) * conj(sig_npk) * sig_nmk;

            fftw_execute_dft_c2r(even_plan, lag_ptr,
                                 W.data() + 2 * n * W.spacing());

            // Odd columns
            lag = state.lambda[m] / (2 * N) * sig_npk * conj(sig_nmk1);

            // TODO Maybe Toole's trick can be applied here. Check that
            fftw_execute(odd_plan);

            // TODO Benchmark with blaze::for_each
            for (int m = 0; m < N; ++m)
                lag[m] *= std::exp(std::complex<double>(0, 1) * (M_PI * m / N));

            column(W_00, 2 * n + 1) = real(lag);
        }

        // W(k,n+N) = (-1)^k W(k,n)
        auto W_01 = submatrix(W, 0, N, N, N);
        for (int n = 0; n < N; ++n) {
            int phase = 1;
            for (int k = 0; k < N; ++k) {
                W_01(k, n) = phase * W_00(k, n);
                phase *= -1;
            }
        }

        // W(k+N,n) = (-1)^n W(k,n)
        auto W_10 = submatrix(W, N, 0, N, N);
        for (int n = 0; 2 * n < N; ++n) {
            column(W_10, 2 * n) = column(W_00, 2 * n);
            column(W_10, 2 * n + 1) = -1 * column(W_00, 2 * n + 1);
        }

        // W(k+N,n+N) = (-1)^(k+n+N) W(k,n)
        int phase = (N % 2 == 0) ? 1 : -1;
        auto W_11 = submatrix(W, N, N, N, N);
        for (int n = 0; n < N; ++n) {
            for (int k = 0; k < N; ++k) {
                W(k + N, n + N) = W(k, n) * phase;
                phase *= -1;
            }
            phase *= -1;
        }
    }

    fftw_destroy_plan(odd_plan);
    fftw_destroy_plan(even_plan);

    return W;
}

}  // namespace Observable

