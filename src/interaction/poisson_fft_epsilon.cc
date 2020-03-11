#include "interaction/poisson_fft_epsilon.h"
#include "parameters.h"
#include "state.h"

#include <cassert>

namespace Poisson {
FFTEpsilon::FFTEpsilon(const Parameters &p, const SimState &state)
    : box(p),
      epsilon{p["Simulation"]["interaction"]["epsilon"]},
      fwd(nullptr),
      bwd(nullptr) {
    // By default, we assume in-place transforms
    auto in_c = reinterpret_cast<const double *>(state.V.data());
    auto out_c = reinterpret_cast<const fftw_complex *>(state.V.data());
    fwd = make_fftw_plan_dft_r2c(box.N, in_c, out_c, FFTW_MEASURE);
    bwd = make_fftw_plan_dft_c2r(box.N, out_c, in_c, FFTW_MEASURE);

    real_ptr = const_cast<double *>(in_c);

    // Nothing to do. For 2D, 3D add fully populated k-grids + rotation
}

void FFTEpsilon::solve(SimState &state) {
    // Compute source of Poisson equation
    state.V = delta_from(state);
    // In place computation, i.e. (i) and (ii) discussed below are FALSE
    solve(state.V, state.V);
}

void FFTEpsilon::solve(blaze::DynamicVector<double> &V,
                       const blaze::DynamicVector<double> &source) {
    // Compute Greens kernel in k-space (no memory allocation!)
    // Note that both odd and even N share the same positive k values
    // if we demand a common upper k-interval boundary. This is different from
    // what is stored in box.kmax. As explained in domain.cc we can shift around
    // as we like because if N is odd the k=-N/2 and k=N/2 mode are redundant
    // (truncating division). We opt for a common upper boundary here to treat N
    // odd and even without additional logic.
    auto kx = blaze::linspace(box.N / 2 + 1, 0.0, box.dk * (box.N / 2));
    auto Ghalf = -1.0 / (kx * kx + epsilon * epsilon);

    // In k-space the complex data is structured as
    //        [ Re(k0), Im(k0), Re(k1), Im(k1), Re(k2), Im(k2), ... ]
    // That is why we have to repeat all values of Ghalf once.
    auto G = kron(Ghalf, blaze::uniform(2, 1.0));

    // Since we perform the FFT without additional memory, we need padding
    // in the target vector. If s.data() == V.data() the procedure is truely
    // in place.
    V.resize(2 * (box.N / 2 + 1));

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
        auto fwd_new =
            make_fftw_plan_dft_r2c(box.N, s_ptr, s_k_ptr, FFTW_MEASURE);
        fftw_execute(fwd_new.get());
    } else {
        fftw_execute_dft_r2c(fwd.get(), s_ptr, s_k_ptr);
    }

    // Solve and normalize in k-space
    V *= 1.0 / box.N * G;

    // Backward FFT is always an in-place transform, so only (ii) needs to be
    // checked
    if (fftw_alignment_of(V_ptr) != fftw_alignment_of(real_ptr)) {
        auto bwd_new =
            make_fftw_plan_dft_c2r(box.N, s_k_ptr, V_ptr, FFTW_MEASURE);
        fftw_execute(bwd_new.get());
    } else {
        fftw_execute_dft_c2r(bwd.get(), s_k_ptr, V_ptr);
    }
    // Strip of padding
    V.resize(box.N);
}

}  // namespace Poisson
