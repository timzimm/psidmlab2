#include "state.h"
#include "cosmology.h"
#include "parameters.h"

SimState::SimState(const Parameters& p)
    : n(0),
      tau(0),
      representation{Representation::Position},
      position_to_momentum{nullptr},
      momentum_to_position{nullptr},
      N_plan{0},
      psi_ptr{nullptr} {}

void SimState::transform(const SimState::Representation target) {
    // Identity transform
    if (target == representation) return;

    auto in = reinterpret_cast<fftw_complex*>(psi.data());
    auto in_real = reinterpret_cast<double*>(in);

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
            make_fftw_plan_dft(N_plan, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
        momentum_to_position =
            make_fftw_plan_dft(N_plan, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    // In particular if we have NOT (i) and NOT (ii) (because the raw arrays are
    // the same or because the raw arrays are different but have the same
    // alignment ) we can reuse the preexisting plan with the new-array
    // interface
    if (target == Representation::Momentum) {
        fftw_execute_dft(position_to_momentum.get(), in, in);
        // FFTW transforms are denormalized
        psi /= N_plan;
    } else {
        fftw_execute_dft(momentum_to_position.get(), in, in);
    }

    representation = target;
};

void operator>>(const SimState& state, Parameters& p) {
    p["Simulation"]["N"] = state.psi.size();
}

