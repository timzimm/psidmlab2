#include "interaction/poisson_fft.h"
#include "parameters.h"
#include "state.h"

namespace Poisson {
FFT::FFT(const Parameters& p)
    : N{p["Simulation"]["N"].get<size_t>()},
      L{p["Simulation"]["L"].get<double>()},
      fft(N / 2 + 1),
      kernel(N / 2 + 1) {
    RCV source_dummy(N);
    forwards.reset(fftw_plan_dft_r2c_1d(
        N, source_dummy.data(), reinterpret_cast<fftw_complex*>(fft.data()),
        FFTW_ESTIMATE));
    backwards.reset(
        fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(fft.data()),
                             source_dummy.data(), FFTW_ESTIMATE));

    // inverse square wavelength and normalization
    for (int k = 1; k < N / 2 + 1; ++k)
        kernel[k] = -1.0 / (4 * M_PI * M_PI * N * k * k / (L * L));

    // makes the DC constraint manifest
    kernel[0] = 0;
}

void FFT::solve(SimState& state) { solve(state.V, delta_from(state)); }

void FFT::solve(blaze::DynamicVector<double>& V,
                const blaze::DynamicVector<double>& s) {
    auto fft_p = reinterpret_cast<fftw_complex*>(fft.data());
    auto source_p =
        const_cast<double*>(reinterpret_cast<const double*>(s.data()));

    fftw_execute_dft_r2c(forwards.get(), source_p, fft_p);

    // Compute fourier coefficients of the potential
    fft = kernel * fft;

    // Get potential by applying the IFFT
    fft_p = reinterpret_cast<fftw_complex*>(fft.data());
    fftw_execute_dft_c2r(backwards.get(), fft_p, V.data());
}
}  // namespace Poisson
