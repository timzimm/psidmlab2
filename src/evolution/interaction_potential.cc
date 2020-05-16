#include "evolution/interaction_potential.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotential::InteractionPotential(const Parameters& p,
                                           const SimState& state,
                                           const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo(cosmo_),
      box(p),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)} {
}

// Non linear phase method with phi_max = pi/2
double InteractionPotential::next_dt(const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau_aux) * max(abs(state.V)));
}

void InteractionPotential::step(SimState& state, const double dt) const {
    state.transform(SimState::Representation::Position);
    const double a = cosmo.a_of_tau(state.tau_aux);
    pot->solve(state);
    state.psi *= exp(-1.0i * a * state.V * dt);

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
