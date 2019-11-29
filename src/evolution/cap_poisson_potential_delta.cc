#include "evolution/cap_poisson_potential_delta.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

CAPPoissonPotentialDelta::CAPPoissonPotentialDelta(const Parameters& p,
                                                   const SimState& state,
                                                   const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      pot{Interaction::make(p["Simulation"]["potential"].get<std::string>(),
                            p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      dx{L / N},
      strength{p["Simulation"]["stepper"]["strength"].get<double>()},
      dt_last{-1},
      CAP(N, strength),
      attenuator(N) {}

// Non linear phase method with phi_max = pi/2
double CAPPoissonPotentialDelta::next_dt(const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void CAPPoissonPotentialDelta::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Position);
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);

    if (dt_last - dt != 0) {
        // 3rd order Taylor expansion in dt of (1-exp(-2CAP*dt)/(2CAP) to avoid
        // cancellation of significant digits in the numerator
        attenuator = dt - dt * dt * CAP + 2.0 / 3 * dt * dt * dt * CAP * CAP;
    }
    pot->solve(state.V, attenuator * real(conj(state.psi) * state.psi) - dt);

    auto drift = exp(-1.0i * a_half * state.V - CAP * dt);
    state.psi *= drift;

    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
