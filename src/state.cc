#include "state.h"
#include "blaze/math/functors/L2Norm.h"
#include "config.h"
#include "cosmology.h"
#include "fftw3.h"
#include "parameters.h"

SimState::SimState(const Domain &box)
    : n(0u), tau(0), dtau(0), tau_aux(0), V(0ul), psi(0ul),
      representation(Representation::Position), position_to_momentum(nullptr),
      momentum_to_position(nullptr), N_plan(0), N_transform(0), norm(0),
      psi_ptr(nullptr) {
    psi.reserve(box.N);
    V.reserve(2 * (box.N / 2 + 1));
    psi.resize(box.N);
    V.resize(box.N);
}

void SimState::transform(const SimState::Representation target) {
    // Identity transform
    if (target == representation)
        return;

    auto in = reinterpret_cast<fftw_complex *>(psi.data());
    auto in_real = reinterpret_cast<double *>(in);

    // Update transformation plans if necessary, i.e.
    //
    // (i) ICs with a different number of spatial points were generated OR
    // (ii) raw pointer changed AND alignment changed
    //
    if (psi.size() != N_plan ||
        (psi_ptr != in_real &&
         fftw_alignment_of(psi_ptr) != fftw_alignment_of(in_real))) {
        N_plan = psi.size();
        psi_ptr = in_real;

        position_to_momentum =
            make_fftw_plan_dft(N_plan, in, in, FFTW_FORWARD, FFTW_MEASURE);
        momentum_to_position =
            make_fftw_plan_dft(N_plan, in, in, FFTW_BACKWARD, FFTW_MEASURE);

        // Cache the norm of the state for which new plans are constructed
        norm = std::real(blaze::l2Norm(psi));
    }

    // In particular if we have NOT (i) and NOT (ii) (because the raw arrays are
    // the same or because the raw arrays are different but have the same
    // alignment ) we can reuse the preexisting plan with the new-array
    // interface
    if (target == Representation::Momentum) {
        N_transform++;
        fftw_execute_dft(position_to_momentum.get(), in, in);
        // FFTW transforms are denormalized
        psi /= N_plan;
    } else {
        N_transform++;
        fftw_execute_dft(momentum_to_position.get(), in, in);
        // Renormalize state to fight numerical errors that crop up after many
        // subsequent FFTs.
        if (N_transform > PSIDMLAB_RENORMALIZATION_THRESHOLD) {
            N_transform = 0;
            psi *= norm / std::real(blaze::l2Norm(psi));
        }
    }

    representation = target;
};
