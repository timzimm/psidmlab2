#include <cassert>
#include <complex>
#include <memory>
#include "poisson_solver.h"

namespace Solvers {
FFT::FFT(const Parameters& p)
    : N(p.N), L(p.L), potential(N), fft(N / 2 + 1), inv_k_sq(N / 2 + 1) {
    RRV source_dummy(N), potential_dummy(N);
    forwards = fftw_plan_dft_r2c_1d(N, source_dummy.data(),
                                    reinterpret_cast<fftw_complex*>(fft.data()),
                                    FFTW_ESTIMATE);
    backwards =
        fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(fft.data()),
                             potential_dummy.data(), FFTW_ESTIMATE);

    auto diag = diagonal(inv_k_sq);
    std::iota(diag.begin() + 1, diag.end(), 1);
    diag = -L * L / (4 * M_PI * M_PI * N) * pow(diag, -2);
    // makes the DC constraint manifest
    diag[0] = 0;
}

FFT::~FFT() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}

FFT::RRV FFT::solve(RRV& source) {
    RRV potential(N);
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
                         potential.data());

    // RVO not applicable. std::move explicetly to avoid deep copy
    return std::move(potential);
}
}  // namespace Solvers
