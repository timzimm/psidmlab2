#include "state.h"

SimState::SimState(const Domain& d)
    : box{d},
      n{0},
      tau{0},
      V(box.N_total),
      psi(box.N_total),
      representation{Representation::Position},
      position_to_momentum{nullptr},
      momentum_to_position{nullptr} {
    auto in = reinterpret_cast<fftw_complex*>(psi.data());

    position_to_momentum.reset(fftw_plan_dft(box.Ns.size(), box.Ns.data(), in,
                                             in, FFTW_FORWARD, FFTW_ESTIMATE));
    momentum_to_position.reset(fftw_plan_dft(box.Ns.size(), box.Ns.data(), in,
                                             in, FFTW_BACKWARD, FFTW_ESTIMATE));
}

void SimState::transform(const SimState::Representation target) {
    // Identity transform
    if (target == representation) return;

    // Due to move semantics and other optimizations, the data adress of blaze's
    // matrices might change in the process. That means that adresses at plan
    // construction time and transformation time might differ. Thus, we
    // explicitly get a pointer to the current raw array and use FFTW advanced
    // interface to perform the transformation
    auto in = reinterpret_cast<fftw_complex*>(psi.data());

    if (target == Representation::Momentum) {
        fftw_execute_dft(position_to_momentum.get(), in, in);
        // FFTW transforms are unnormalized
        psi /= box.N_total;
    } else {
        fftw_execute_dft(momentum_to_position.get(), in, in);
    }

    representation = target;
};
