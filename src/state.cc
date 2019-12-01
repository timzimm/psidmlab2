#include "state.h"
#include "cosmology.h"
#include "logging.h"
#include "parameters.h"

#include <functional>
#include <numeric>

SimState::SimState(const Parameters& p)
    : n{0},
      tau{0},
      dofs{p["Domain Properties"]["dof"].get<std::vector<int>>()},
      box{},
      dims{dofs.size()},
      N_total{0},
      representation{Representation::Position},
      position_to_momentum{nullptr},
      momentum_to_position{nullptr} {
    for (auto& limit : p["Domain Properties"]["limits"]) {
        box.emplace_back(limit.get<interval_t>());
    }
    if (box.size() != dims) {
        std::cout << ERRORTAG("Ill-defined box") << std::endl;
        exit(1);
    }

    N_total = std::accumulate(dofs.begin(), dofs.end(), 1, std::multiplies<>());
    psi.resize(N_total);
    V.resize(N_total);

    auto in = reinterpret_cast<fftw_complex*>(psi.data());
    position_to_momentum = fftw_plan_dft(dofs.size(), dofs.data(), in, in,
                                         FFTW_FORWARD, FFTW_ESTIMATE);
    momentum_to_position = fftw_plan_dft(dofs.size(), dofs.data(), in, in,
                                         FFTW_BACKWARD, FFTW_ESTIMATE);
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
        fftw_execute_dft(position_to_momentum, in, in);
        // FFTW transforms are denormalized
        psi /= N_total;
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
    p["Simulation"]["N"] = state.psi.size();
}

