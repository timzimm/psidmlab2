#include "evolution/interaction_potential_magnus.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotentialMagnus::InteractionPotentialMagnus(const Parameters& p,
                                                       const SimState& state,
                                                       const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo(cosmo_),
      box(p),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)} {
}

// Non linear phase method
double InteractionPotentialMagnus::next_dt(const SimState& state) const {
    const double phi_max = 0.1;
    return phi_max / (cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void InteractionPotentialMagnus::step(SimState& state, const double dt) const {
    state.transform(SimState::Representation::Position);
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    pot->solve(state);
    state.psi *= exp(-1.0i * a_half * state.V * dt);

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
