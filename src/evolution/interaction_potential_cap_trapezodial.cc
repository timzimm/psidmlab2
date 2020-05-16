#include "evolution/interaction_potential_cap_trapezodial.h"
#include "blaze/math/views/Subvector.h"
#include "parameters.h"

#include "cosmology.h"

#include <string>

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotentialCAPTrapezodial::InteractionPotentialCAPTrapezodial(
    const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      box{p},
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)},
      dt_last{-1},
      CAP(box.N, 1),
      attenuator(box.N) {
    const double strength = p["Simulation"]["stepper"]["strength"];
    const double width = p["Simulation"]["stepper"]["width"];
    const double rp = p["Simulation"]["stepper"]["rp"];
    // In leading order tanh(1.0/width x) = 1.0/width * x
    // So dx = dy / dy/dx = 5 * width (5 because tanh is bounded in [-1, 1] and
    // more than doubling that should make the transition smooth enough)
    const double rs = rp + 5 * width;

    auto x = linspace(box.N, box.xmin, box.xmax);
    CAP *= strength / 2.0 * (tanh(1.0 / width * (abs(x) - rs)) + 1);
    // Heaviside function. CAP is zero inside the pysical domain.
    int rp_N = static_cast<int>(rp / box.dx + 0.5);
    subvector(CAP, box.N / 2 - rp_N, 2 * rp_N) = 0;
}

// Non linear phase method with phi_max = pi/2
double InteractionPotentialCAPTrapezodial::next_dt(
    const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void InteractionPotentialCAPTrapezodial::step(SimState& state,
                                              const double dt) const {
    state.transform(SimState::Representation::Position);

    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    if (dt_last - dt != 0) {
        attenuator = exp(-CAP * dt / 2.0);
    }

    // In place computation of potential
    state.V = attenuator * rho_from(state);
    pot->solve(state.V, state.V);

    state.psi *= exp(-1.0i * a_half * state.V * dt) * attenuator;

    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
