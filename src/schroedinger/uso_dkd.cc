#include <cassert>

#include "cosmology.h"
#include "parameters.h"
#include "schroedinger/uso_dkd.h"

namespace Schroedinger {
using namespace blaze;
using namespace std::complex_literals;

USO_DKD::USO_DKD(const Parameters& p, const SimState& state,
                 const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      dt_last{-1},
      a{cosmo.a_of_tau(-1)},
      k_squared(N),
      kick(N, state.M) {
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.end() - (N + 1) / 2,
                k_squared.end());
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;
}

double USO_DKD::next_dt(const SimState& state) const {
    return std::min(
        {L * L / (2 * N * N * M_PI), M_PI / (2 * a * max(abs(state.V)))});
}

void USO_DKD::step(SimState& state, const double dt) {
    CCM& psis = state.psis;
    const double a_next = cosmo.a_of_tau(state.tau + dt);

    // Drift step
    state.transform(SimState::Representation::Position);
    psis =
        expand(exp(-0.5i * 0.5 * (a + a_next) * state.V * dt), state.M) % psis;

    state.transform(SimState::Representation::Momentum);

    if (dt_last - dt != 0) {
        kick = expand(exp(-0.5i * k_squared * dt), state.M);
    }
    psis = kick % psis;

    state.transform(SimState::Representation::Position);
    pot->solve(state);

    psis =
        expand(exp(-0.5i * 0.5 * (a + a_next) * state.V * dt), state.M) % psis;

    // state is now @ tau + dt
    pot->solve(state);
    state.tau += dt;
    a = a_next;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
