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
      M_plan{0} {}

void SimState::transform(const SimState::Representation target) {
    // Idenity transform
    if (target == representation) return;

    // Due to move semantics and other optimizations, the data adress of blaze's
    // matrices might change in the process. That means that adresses at plan
    // construction time and transformation time might differ. Thus, we
    // explicitly get a pointer to the current raw array and use FFTW advanced
    // interface to perform the transformation
    auto in = reinterpret_cast<fftw_complex*>(psis.data());

    // Update transformation plans if necessary, e.g. if ICs with a different
    // number of wavefunctions or spatial points was generated
    if (psis.columns() != M_plan || psis.rows() != N_plan) {
        M = M_plan = psis.columns();
        N_plan = psis.rows();
        fftw_destroy_plan(momentum_to_position);
        int offset = psis.spacing();

        position_to_momentum =
            fftw_plan_many_dft(1, &N_plan, M_plan, in, nullptr, 1, offset, in,
                               nullptr, 1, offset, FFTW_FORWARD, FFTW_ESTIMATE);
        momentum_to_position = fftw_plan_many_dft(
            1, &N_plan, M_plan, in, nullptr, 1, offset, in, nullptr, 1, offset,
            FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    if (target == Representation::Momentum) {
        fftw_execute_dft(position_to_momentum, in, in);
        // FFTW transforms are denormalized
        psis /= N_plan;
    } else {
        fftw_execute_dft(momentum_to_position, in, in);
    }

    representation = target;
};

SimState::~SimState() {
    fftw_destroy_plan(position_to_momentum);
    fftw_destroy_plan(momentum_to_position);
}

void operator>>(const SimState& state, Parameters& p) {
    p["Simulation"]["N"] = state.psis.rows();
    p["Initial Conditions"]["M"] = state.psis.columns();
}

