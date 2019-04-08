#include "poisson/fft.h"
#include <cassert>
#include "common.h"

namespace Poisson {
FFT::FFT(const Parameters& p)
    : N(p.N), L(p.L), fft(N / 2 + 1), inv_k_sq(N / 2 + 1) {
    RCV source_dummy(N);
    RCV potential_dummy(N);
    forwards = fftw_plan_dft_r2c_1d(N, source_dummy.data(),
                                    reinterpret_cast<fftw_complex*>(fft.data()),
                                    FFTW_ESTIMATE);
    backwards =
        fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(fft.data()),
                             potential_dummy.data(), FFTW_ESTIMATE);

    auto diag = blaze::diagonal(inv_k_sq);
    for (int k = 1; k < N / 2 + 1; ++k)
        diag[k] = -L * L / (4 * M_PI * M_PI * N) / (k * k);
    // makes the DC constraint manifest
    diag[0] = 0;
}

FFT::~FFT() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}

void FFT::operator()(SimState& state) {
    // Calculate source term
    auto psi2 = blaze::real(state.psis % state.psis);
    RCV source = psi2 * state.lambda;

    // fftw_execute_dft_r2c and its inverse are applicable in this situation
    // because blaze takes care of proper alignment
    fftw_execute_dft_r2c(forwards, source.data(),
                         reinterpret_cast<fftw_complex*>(fft.data()));

    // Check if the DC constrained is satisfied
    assert(std::abs(fft[0]) < epsilon);

    // Compute fourier coefficients of the potential
    fft = inv_k_sq * fft;

    // Get potential by applying the IFFT
    fftw_execute_dft_c2r(backwards, reinterpret_cast<fftw_complex*>(fft.data()),
                         state.V.data());
}
}  // namespace Poisson
