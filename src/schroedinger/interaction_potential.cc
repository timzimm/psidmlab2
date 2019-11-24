#include "schroedinger/interaction_potential.h"
#include "cosmology.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotential::InteractionPotential(const Parameters& p,
                                           const SimState& state,
                                           const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)} {}

// Non linear phase method with phi_max = pi/2
double InteractionPotential::next_dt(const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void InteractionPotential::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Position);
    // |psi|^2 is invariant. Thus we can use the input state to determine the
    // potential for the entire time step.
    pot->solve(state);
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    auto drift = expand(exp(-1.0i * a_half * state.V * dt), state.M);
    state.psis = drift % state.psis;

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
