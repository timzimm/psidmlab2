#include "evolution/kinetic.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>

namespace Schroedinger {
using namespace blaze;
using namespace std::complex_literals;

Kinetic::Kinetic(const Parameters& p, const SimState& state, const Cosmology&)
    : DefaultDriver(p),
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      dt_last{-1},
      kx2(N) {
    // Next largest even number
    const int NN = (N % 2) ? N + 1 : N;
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    //          f0 f1 ... fN/2, f-N/2+1 ... f-1
    // This rotation makes it hard to generate a blaze expression at compile
    // time. If we still want to get around a raw loop it seems best to do an
    // allocation for the rotated dimension. In higher dimensions only one axis
    // is rotated meaning we can use linspace for all other axes.
    std::iota(kx2.begin(), kx2.end(), -NN / 2 + 1);
    std::rotate(kx2.begin(), kx2.begin() + NN / 2 - 1, kx2.end());
    kx2 = 4 * M_PI * M_PI / (L * L) * kx2 * kx2;
}

// Limit max phase change to pi/2 per step
double Kinetic::next_dt(const SimState& state) const {
    return L * L / (N * N * M_PI);
}

void Kinetic::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Momentum);

    // Since in 1D we store the k values anyways,
    // we only need to update the matrix exponential if an altered time step
    // size is required
    if (dt - dt_last != 0) {
        state.psi *= exp(-0.5i * kx2 * dt);
    }
    // in d > 1 we "recompute" the k-grid in each step to get around memory
    // allocation. Note that using a raw loop does exactly the same but less
    // expressive.

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
