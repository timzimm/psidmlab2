#include "schroedinger/kinetic.h"
#include "parameters.h"

#include <algorithm>

namespace Schroedinger {
using namespace blaze;
using namespace std::complex_literals;

Kinetic::Kinetic(const Parameters& p, const SimState& state, const Cosmology&)
    : DefaultDriver(p),
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      dt_last{-1},
      k_squared(N),
      kick(N, state.M) {
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.end() - (N + 1) / 2,
                k_squared.end());
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;
}

// Limit max phase change to pi/2 per step
double Kinetic::next_dt(const SimState& state) const {
    return L * L / (N * N * M_PI);
}

void Kinetic::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Momentum);

    // We only need to update the matrix exponential if an altered time step
    // size is required
    if (dt - dt_last != 0) {
        kick = expand(exp(-0.5i * k_squared * dt), state.M);
    }

    state.psis = kick % state.psis;

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
