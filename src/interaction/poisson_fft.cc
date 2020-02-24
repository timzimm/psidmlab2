#include "interaction/poisson_fft.h"
#include "fftw.h"
#include "fftw3.h"
#include "parameters.h"
#include "state.h"

namespace Poisson {
FFT::FFT(const Parameters &p, const SimState &state)
    : N{p["Simulation"]["N"].get<size_t>()},
      L{p["Simulation"]["L"].get<double>()} {
    // By default, we assume in-place transforms
    auto in_c = reinterpret_cast<const double *>(state.V.data());
    auto out_c = reinterpret_cast<const fftw_complex *>(state.V.data());

    // FFTW does not know about constness
    real_ptr = const_cast<double *>(in_c);
    auto out = const_cast<fftw_complex *>(out_c);

    fftw_plan_ptr fwd(fftw_plan_dft_r2c_1d(N, real_ptr, out, FFTW_ESTIMATE));
    fftw_plan_ptr bwd(fftw_plan_dft_c2r_1d(N, out, real_ptr, FFTW_ESTIMATE));

    // Nothing to do. For 2D, 3D add fully populated k-grids + rotation
}

void FFT::solve(SimState &state) {
    // Compute source of Poisson equation
    state.V = delta_from(state);
    // In place computation, i.e. (i) and (ii) discussed below are FALSE
    solve(state.V, state.V);
}

void FFT::solve(blaze::DynamicVector<double> &V,
                const blaze::DynamicVector<double> &source) {
    // Compute Greens kernel in k-space (no memory allocation!)
    // We omit k = 0
    auto kx = blaze::linspace(N / 2, 2 * M_PI / L, M_PI / (L / N));
    auto Ghalf = -1.0 / (kx * kx);
    // In k-space the complex data is structured as
    //        [ Re(k1), Im(k1), Re(k2), Im(k2), ... ]
    // That is why we have to repeat all values of Ghalf once.
    auto G = kron(Ghalf, blaze::uniform(2, 1.0));

    // Since we perform the FFT without additional memory, we need padding in
    // the target vector. If source.data() == V.data() the procedure is truely
    // in place.
    V.resize(2 * (N / 2 + 1));
    // FFTW doesn't know about constness
    auto s_ptr = const_cast<double *>(source.data());
    auto s_k_ptr = reinterpret_cast<fftw_complex *>(V.data());
    auto V_ptr = V.data();

    // We only need to replan the forward FFT if
    // (i) V.data() != source.data() (transform is out-of-place) OR
    // (ii) V.data() has different alignment than state.V
    // This saves us the planing time which might be significant for small
    // problems
    if (s_ptr != V_ptr ||
        fftw_alignment_of(s_ptr) != fftw_alignment_of(real_ptr)) {
        fftw_plan_ptr fwd_new(
            fftw_plan_dft_r2c_1d(N, s_ptr, s_k_ptr, FFTW_ESTIMATE));
        fftw_execute(fwd_new.get());
    } else {
        // New array execute functions. See:
        // http://www.fftw.org/fftw3_doc/New_002darray-Execute-Functions.html
        fftw_execute_dft_r2c(fwd.get(), s_ptr, s_k_ptr);
    }

    // Solve and normalize in k-space
    subvector(V, 2, V.size() - 2) *= 1.0 / N * G;
    // Assume DC mode vanishes (even if it does not!) and singularity is
    // irrelevenat
    V[0] = 0;
    V[1] = 0;

    // Backward FFT is always an in-place transform, so only (ii) needs to be
    // checked
    if (fftw_alignment_of(V_ptr) != fftw_alignment_of(real_ptr)) {
        fftw_plan_ptr bwd_new(
            fftw_plan_dft_c2r_1d(N, s_k_ptr, V_ptr, FFTW_ESTIMATE));
        fftw_execute(bwd_new.get());
    } else {
        fftw_execute_dft_c2r(bwd.get(), s_k_ptr, V_ptr);
    }
    // Strip of padding
    V.resize(N);
}

}  // namespace Poisson
