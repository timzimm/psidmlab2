#include "schroedinger/cap_poisson_potential.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

CAPPoissonPotential::CAPPoissonPotential(const Parameters& p,
                                         const SimState& state,
                                         const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      dx{L / N},
      strength{p["Simulation"]["stepper"]["strength"].get<double>()},
      width{p["Simulation"]["stepper"]["width"].get<double>()},
      dt_last{-1},
      CAP(N),
      attenuator(N) {
    std::iota(CAP.begin(), CAP.end(), -N / 2);
    CAP = strength * (pow(cosh(1.0 / width * (CAP * dx - L / 2)), -2) +
                      pow(cosh(1.0 / width * (CAP * dx + L / 2)), -2));
}

double CAPPoissonPotential::next_dt(const SimState& state) const {
    return M_PI / (cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void CAPPoissonPotential::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Position);

    if (dt_last - dt != 0) {
        // 3rd order Taylor expansion in dt of (1-exp(-2CAP*dt)/(2CAP) to avoid
        // cancellation of significant digits in the numerator
        attenuator = dt - dt * dt * CAP + 2.0 / 3 * dt * dt * dt * CAP * CAP;
    }
    pot->solve(
        state.V,
        attenuator * (real(conj(state.psis) % state.psis) * state.lambda) - dt);

    const double a = cosmo.a_of_tau(state.tau);
    const double a_next = cosmo.a_of_tau(state.tau + dt);
    auto drift =
        expand(exp(-1.0i * 0.5 * (a + a_next) * state.V - CAP * dt), state.M);
    state.psis = drift % state.psis;

    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
