#include "state.h"
#include "blaze/math/functors/L2Norm.h"
#include "config.h"
#include "cosmology.h"
#include "fftw3.h"
#include "logging.h"
#include "parameters.h"

SimState::SimState(const Domain& box_)
    : n(0u),
      tau(0),
      tau_aux(0),
      V(0ul),
      psi(0ul),
      box(box_),
      representation(Representation::Position),
      position_to_momentum(nullptr),
      momentum_to_position(nullptr),
      N_plan(box.N),
      N_transform(0),
      psi_ptr(nullptr) {
    psi.reserve(box.N);
    psi.resize(box.N);
    // This allows for in place FFT transforms in state.V without the need of
    // reallocations
    V.reserve(2 * (box.N / 2 + 1));
    V.resize(box.N);

    auto in = reinterpret_cast<fftw_complex*>(psi.data());
    psi_ptr = reinterpret_cast<double*>(in);

    if (box.bc == Domain::BoundaryCondition::Periodic) {
        position_to_momentum =
            make_fftw_plan_dft(N_plan, in, in, FFTW_FORWARD, FFTW_MEASURE);
        momentum_to_position =
            make_fftw_plan_dft(N_plan, in, in, FFTW_BACKWARD, FFTW_MEASURE);
    } else if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
        // Under Homogeneous Dirichlet Conditions, we transform real and
        // imaginary part by a DST-1. This is achieved by directly accessing
        // the raw array of the complex psi vector and acting on it with two
        // independent DSTs of stride 2.
        const int howmany = 2;
        const int n[] = {N_plan, N_plan};
        const int rank = 1;
        const int idist = 1;
        const int odist = idist;
        const int istride = 2;
        const int ostride = istride;
        const fftw_r2r_kind kinds[] = {FFTW_RODFT00, FFTW_RODFT00};
        position_to_momentum = make_fftw_plan_many_r2r(
            rank, n, howmany, psi_ptr, nullptr, istride, idist, psi_ptr,
            nullptr, ostride, odist, kinds, FFTW_MEASURE);
        momentum_to_position = make_fftw_plan_many_r2r(
            rank, n, howmany, psi_ptr, nullptr, istride, idist, psi_ptr,
            nullptr, ostride, odist, kinds, FFTW_MEASURE);
    }
}

void SimState::transform(const SimState::Representation target) {
    static double norm = std::real(blaze::l2Norm(psi));
    // Identity transform
    if (target == representation) return;

#ifndef NDEBUG
    auto in = reinterpret_cast<fftw_complex*>(psi.data());
    auto in_real = reinterpret_cast<double*>(in);
    // Check if raw pointer changed AND alignment changed
    if (psi.size() != N_plan ||
        (psi_ptr != in_real &&
         fftw_alignment_of(psi_ptr) != fftw_alignment_of(in_real))) {
        std::cout << ERRORTAG("Plan missmatch") << std::endl;
        exit(1);
    }
#endif

    // In particular if we have NOT (i) and NOT (ii) (because the raw arrays are
    // the same or because the raw arrays are different but have the same
    // alignment ) we can reuse the preexisting plan with the new-array
    // interface
    if (target == Representation::Momentum) {
        N_transform++;
        /* fftw_execute_dft(position_to_momentum.get(), in, in); */
        fftw_execute(position_to_momentum.get());
        // FFTW transforms are denormalized
        if (box.bc == Domain::BoundaryCondition::Periodic) {
            psi /= N_plan;
        } else if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
            psi /= 2 * (N_plan + 1.0);  // see FFTW manual
        }
    } else {
        N_transform++;
        /* fftw_execute_dft(momentum_to_position.get(), in, in); */
        fftw_execute(momentum_to_position.get());
        // Renormalize state to fight numerical errors that crop up after many
        // subsequent FFTs.
        if (N_transform > PSIDMLAB_RENORMALIZATION_THRESHOLD) {
            N_transform = 0;
            psi *= norm / std::real(blaze::l2Norm(psi));
        }
    }

    representation = target;
};
