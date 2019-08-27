#include "schroedinger/poisson_potential.h"
#include "cosmology.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

PoissonPotential::PoissonPotential(const Parameters& p, const SimState& state,
                                   const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)} {}

// Non linear phase method with phi_max = pi/2
double PoissonPotential::next_dt(const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void PoissonPotential::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Position);
    // |psi|^2 is invariant. Thus we can use the input state to determine the
    // potential for the entire time step.
    pot->solve(state);
    const double a = cosmo.a_of_tau(state.tau);
    const double a_next = cosmo.a_of_tau(state.tau + dt);
    auto drift =
        expand(exp(-1.0i * 0.5 * (a + a_next) * state.V * dt), state.M);
    state.psis = drift % state.psis;

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
