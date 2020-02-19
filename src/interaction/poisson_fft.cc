#include "interaction/poisson_fft.h"
#include "fftw.h"
#include "parameters.h"
#include "state.h"

namespace Poisson {
FFT::FFT(const Parameters &p)
    : N{p["Simulation"]["N"].get<size_t>()},
      L{p["Simulation"]["L"].get<double>()} {
    // Nothing to do here. For 2D, 3D add fully populated k-grids + rotation
}

void FFT::solve(SimState &state) {
    // Compute source of Poisson equation
    state.V = delta_from(state);
    // In place computation
    solve(state.V, state.V);
}

void FFT::solve(blaze::DynamicVector<double> &V,
                const blaze::DynamicVector<double> &source) {
    // Compute Greens kernel in k-space (no memory allocation!)
    auto kx = blaze::linspace(N / 2 + 1, 0.0, M_PI / (L / N));
    auto Ghalf = -1.0 / (kx * kx);
    // In k-space the complex data is structured as Re(1) Im(1) Re(2) Im(2)...
    // That is why we have to repeat all values of Ghalf once.
    auto G = kron(Ghalf, blaze::uniform(2, 1.0));

    // Since we perform the FFT without additional memory, we need padding in
    // the target vector. If s.data() == V.data() the procedure is truely in
    // place.
    V.resize(2 * (N / 2 + 1));
    // FFTW doesn't know about constness
    auto s_ptr = const_cast<double *>(source.data());
    auto s_k_ptr = reinterpret_cast<fftw_complex *>(V.data());
    auto V_ptr = V.data();

    // This might be an in-place transform (if s.data() == V.data())
    fftw_plan_ptr fwd(fftw_plan_dft_r2c_1d(N, s_ptr, s_k_ptr, FFTW_ESTIMATE));
    // This is always an in-place transform
    fftw_plan_ptr bwd(fftw_plan_dft_c2r_1d(N, s_k_ptr, V_ptr, FFTW_ESTIMATE));

    fftw_execute(fwd.get());

    // Solve and normalize in k-space
    V *= 1.0 / N * G;
    // Assume DC mode vanishes (even if it does not!) and singularity is
    // irrelevenat
    V[0] = 0;
    V[1] = 0;

    fftw_execute(bwd.get());
    // Strip of padding
    V.resize(N);
}

}  // namespace Poisson
