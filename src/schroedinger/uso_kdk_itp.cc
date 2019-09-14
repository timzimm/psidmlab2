#include "schroedinger/uso_kdk_itp.h"
#include "cosmology.h"
#include "parameters.h"

#include <algorithm>

using namespace blaze;
using namespace std::complex_literals;

namespace Schroedinger {

USO_KDK_ITP::USO_KDK_ITP(const Parameters& p, const SimState& state,
                         const Cosmology& cosmo_)
    : DefaultDriver{p},
      cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      a{cosmo.a_of_tau(-1)},
      k_squared(N),
      kick(N, state.M),
      dt_last{-1} {
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.end() - (N + 1) / 2,
                k_squared.end());
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;
}

void USO_KDK_ITP::step(SimState& state, const double dt) {
    CCM& psis = state.psis;
    const double a_next = cosmo.a_of_tau(state.tau + dt);

    // Set kick operator once for the entire time step
    if (dt - dt_last != 0) {
        kick = expand(exp(-0.5 * 0.5 * k_squared * dt), state.M);
    }

    state.transform(SimState::Representation::Momentum);
    psis = kick % psis;

    state.transform(SimState::Representation::Position);
    // Recompute potential
    pot->solve(state);

    // Set drift operator with intermediate potential
    auto drift = expand(exp(-1.0 * 0.5 * (a + a_next) * state.V * dt), state.M);
    psis = drift % psis;

    state.transform(SimState::Representation::Momentum);
    psis = kick % psis;
    state.transform(SimState::Representation::Position);

    double norm2 =
        L / N * sum(real(conj(state.psis) % state.psis) * state.lambda);

    column(state.psis, 0) /= sqrt(norm2);

    state.tau += dt;
    a = a_next;
    dt_last = dt;
    state.n += 1;
}

double USO_KDK_ITP::next_dt(const SimState& state) const {
    return std::min({L * L / (N * N * M_PI), M_PI / (a * max(abs(state.V)))});
}

}  // namespace Schroedinger
