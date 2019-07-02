#include "poisson/fft.h"
#include "parameters.h"
#include "state.h"

#include <cassert>

namespace Poisson {
FFT::FFT(const Parameters& p)
    : N(p["Simulation"]["N"].get<size_t>()),
      L(p["Simulation"]["L"].get<double>()),
      fft(N / 2 + 1),
      source(N),
      inv_k_sq(N / 2 + 1) {
    RCV potential_dummy(N);
    forwards = fftw_plan_dft_r2c_1d(N, source.data(),
                                    reinterpret_cast<fftw_complex*>(fft.data()),
                                    FFTW_ESTIMATE);
    backwards =
        fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(fft.data()),
                             potential_dummy.data(), FFTW_ESTIMATE);

    auto diag = blaze::diagonal(inv_k_sq);

    // inverse square wavelength and normalization
    for (int k = 1; k < N / 2 + 1; ++k)
        diag[k] = -L * L / (4 * M_PI * M_PI * N) / (k * k);

    // makes the DC constraint manifest
    diag[0] = 0;
}

FFT::~FFT() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}

void FFT::solve(SimState& state) {
    source = delta_from(state);

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
