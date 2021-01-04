#include "interaction/periodic_convolution.h"
#include "fftw3.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

#include <gsl/gsl_sf_hyperg.h>

PeriodicConvolution::PeriodicConvolution(const Parameters &p,
                                         const SimState &state)
    : box(p), G_k(2 * (box.N / 2)), fwd(nullptr), bwd(nullptr) {
    // By default, we assume in-place transforms
    auto in_c = reinterpret_cast<const double *>(state.V.data());
    auto out_c = reinterpret_cast<const fftw_complex *>(state.V.data());

    std::cout << INFOTAG("Planning In-place Interaction FFTs...") << std::flush;
    fwd = make_fftw_plan_dft_r2c(box.N, in_c, out_c, FFTW_PATIENT);
    bwd = make_fftw_plan_dft_c2r(box.N, out_c, in_c, FFTW_PATIENT);
    std::cout << "...done" << std::endl;

    // Compute and store Greens kernel in k-space
    // Note that both odd and even N share the same positive k values
    // if we demand a common upper k-interval boundary. This is different from
    // what is stored in box.kmax. As explained in domain.cc we can shift around
    // as we like because if N is even the k=-N/2 and k=N/2 mode are redundant.
    //
    // We opt for a common upper boundary here to treat N odd and even without
    // additional logic.
    //
    // k-grid construction:
    // (i)  We omit k=0 explicitely to make the DC-constraint mannifest
    // (ii) In k-space the complex data is structured as
    //        [ Re(k1), Im(k1), Re(k2), Im(k2), ... ]
    //      That is why we have to repeat all values of k once.
    auto k = kron(blaze::linspace(box.N / 2, box.dk, box.dk * (box.N / 2)),
                  blaze::uniform(2, 1.0));

    const InteractionType type = p["Simulation"]["interaction"]["type"];
    if (type == InteractionType::Poisson) {
        G_k = -1.0 / (k * k);
    } else if (type == InteractionType::ScreenedPoisson) {
        const double eps = p["Simulation"]["interaction"]["epsilon"];
        G_k = -1.0 / (k * k + eps * eps);
    } else if (type == InteractionType::LAM) {
        const double eps = p["Simulation"]["interaction"]["epsilon"];
        G_k = map(k, [&eps](double k) {
            return -1.0 / (4 * M_PI) *
                   gsl_sf_hyperg_U_int(1, 1, 0.5 * eps * eps * k * k);
        });
    } else {
        std::cerr << ERRORTAG("Interaction not supported") << std::endl;
        exit(1);
    }

    real_ptr = const_cast<double *>(in_c);
}

void PeriodicConvolution::solve(SimState &state) {
    state.V = rho_from(state);
    // In place computation, i.e. (i) and (ii) discussed below are FALSE
    solve(state.V, state.V);
}

void PeriodicConvolution::solve(blaze::DynamicVector<double> &V,
                                const blaze::DynamicVector<double> &source) {
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
    if (s_ptr != V_ptr ||
        fftw_alignment_of(s_ptr) != fftw_alignment_of(real_ptr)) {

        std::cout << INFOTAG("Re-planning Forward Interaction FFT...")
                  << std::flush;
        auto fwd_new =
            make_fftw_plan_dft_r2c(box.N, s_ptr, s_k_ptr, FFTW_ESTIMATE);
        std::cout << "...done" << std::endl;

        fftw_execute(fwd_new.get());
    } else {
        fftw_execute_dft_r2c(fwd.get(), s_ptr, s_k_ptr);
    }

    // Convolution Theorem + Normalization
    subvector(V, 2, V.size() - 2) *= 1.0 / box.N * G_k;

    // Assume DC mode vanishes (even if it does not!)
    V[0] = 0;
    V[1] = 0;

    // Backward FFT is always an in-place transform, so only (ii) needs to be
    // checked
    if (fftw_alignment_of(V_ptr) != fftw_alignment_of(real_ptr)) {

        std::cout << INFOTAG("Re-planning Backward Interaction FFT...")
                  << std::flush;
        auto bwd_new =
            make_fftw_plan_dft_c2r(box.N, s_k_ptr, V_ptr, FFTW_ESTIMATE);
        std::cout << "...done" << std::endl;
        fftw_execute(bwd_new.get());
    } else {
        fftw_execute_dft_c2r(bwd.get(), s_k_ptr, V_ptr);
    }
    // Strip of padding
    V.resize(box.N);
}
