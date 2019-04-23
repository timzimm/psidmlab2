#include "observables.h"
#include <complex>
#include <iostream>
#include "blaze/math/UniformVector.h"
#include "parameters.h"
#include "state.h"

#include <blaze/math/Rows.h>
#include <blaze/math/Submatrix.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
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

auto peyrin(const blaze::DynamicVector<std::complex<double>>& signal) {
    using namespace blaze;
    const int N = signal.size();
    DynamicMatrix<double, columnMajor> W(2 * N, 2 * N);
    // We only compute the first quadrant of W. All other quadrants follow by
    // symmetry properties of W.
    auto W_00 = blaze::submatrix(W, 0, 0, N, N);

    DynamicVector<std::complex<double>> strip(N);

    auto in_out = reinterpret_cast<fftw_complex*>(strip.data());

    fftw_plan one =
        fftw_plan_dft_1d(N, in_out, in_out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int n = 0; 2 * n < N; ++n) {
        // Even
        for (int k = 0; k < N; ++k) {
            int kpn = (k + n) % N;
            int nmk = (n - k) % N;
            nmk = (nmk < 0) ? nmk + N : nmk;
            strip[k] = signal[kpn] * std::conj(signal[nmk]);
        }
        fftw_execute(one);
        column(W_00, 2 * n) = 1.0 / (2 * N) * real(strip);
        // Odd
        for (int k = 0; k < N; ++k) {
            int kpn = (k + n) % N;
            int nmk = (n - k + 1) % N;
            nmk = (nmk < 0) ? nmk + N : nmk;
            strip[k] = signal[kpn] * std::conj(signal[nmk]);
        }
        fftw_execute(one);
        for (int m = 0; m < N; ++m)
            strip[m] *= std::exp(std::complex<double>(0, 1) * (M_PI * m / N));

        column(W_00, 2 * n + 1) = 1.0 / (2 * N) * real(strip);
    }
    auto W_01 = submatrix(W, 0, N, N, N) = W_00;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; 2 * k < N; ++k) {
            W_01(2 * k + 1, n) *= -1;
        }
    }
    auto W_10 = submatrix(W, N, 0, N, N) = W_00;
    for (int n = 0; 2 * n < N; ++n) {
        for (int k = 0; 2 * k < N; ++k) {
            W_10(k, 2 * n + 1) *= -1;
        }
    }
    auto W_11 = submatrix(W, N, N, N, N) = W_00;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            W_11(k, n) *= std::pow(-1, k + n + N);
        }
    }

    return W;
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

    // Load chirp data
    std::ifstream input("/tmp/chirp/chirp.txt");
    std::copy(std::istream_iterator<std::complex<double>>(input),
              std::istream_iterator<std::complex<double>>(), chirp.begin());

    // Phase factor depends on signal size
    auto rho = [this](int tau, int nu) {
        const std::complex<double> i(0, 1);
        int tau_nu = (tau * nu) % N;
        if (N % 2 == 0)
            if (tau_nu < N / 2)
                return std::exp(i * (M_PI / N * tau_nu));
            else
                return std::exp(i * (M_PI / N * (tau_nu - N)));
        int inverse = (N + 1) / 2;
        return std::exp(i * (2 * M_PI / N * inverse * tau_nu));
    };

    // Ambuigity matrix
    blaze::DynamicMatrix<std::complex<double>> A(N, N, 0);

    auto& signal = chirp;

    auto in_out = reinterpret_cast<fftw_complex*>(blaze::row(A, 0).data());
    fftw_plan one =
        fftw_plan_dft_1d(N, in_out, in_out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan two =
        fftw_plan_dft_2d(N, N, in_out, in_out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Construct Ambuigity matrix columwise (modulo the phase factor rho)
    for (int tau = 0; tau < N; ++tau) {
        for (int l = 0; l < N; ++l) {
            A(tau, l) = 1.0 / N * std::conj(signal[(l + tau) % N]) * signal[l];
        }
        fftw_execute_dft(one, in_out, in_out);
        in_out = reinterpret_cast<fftw_complex*>(blaze::row(A, tau).data());
    }
    // At this point we constructed the ambiguity matrix. The Wigner
    // distribution is given by W = IDFT(A). We still have to take rho into
    // account and add an additional phase to move the frequency origin into
    // the matrix center.
    const std::complex<double> i(0, 1);
    for (int tau = 0; tau < N; ++tau) {
        for (int nu = 0; nu < N; ++nu) {
            A(tau, nu) *= std::exp((tau * nu) % N * M_PI / N * i);
        }
    }
    fftw_execute(two);
    blaze::DynamicMatrix<double, blaze::columnMajor> W = blaze::real(A);

    return peyrin(signal);
}  // namespace Observable

}  // namespace Observable

