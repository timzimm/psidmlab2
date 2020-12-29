#include "evolution/kinetic.h"
#include "parameters.h"

#include <algorithm>
#include <cmath>

namespace Schroedinger {
using namespace blaze;
using namespace std::complex_literals;

Kinetic::Kinetic(const Parameters &p, const SimState &state, const Cosmology &)
    : DefaultDriver(p), box(p), dt_last{-1}, k2_max(M_PI / box.dx), kx2(box.N),
      U(box.N) {
    // Next even number NN at least N
    const int NN = (box.N % 2) ? box.N + 1 : box.N;
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    //          f_0 f_1 ... f_box.N/2, f_-NN/2+1 ... f_-1
    // This rotation makes it hard to generate a blaze expression at compile
    // time. If we still want to get around a raw loop it seems best to:
    //
    // (i) do an allocation for the rotated dimension.
    // (ii) have a complex vector for the evolution operator
    //
    // (i) turns (ii) into a performance optimization because we can omit a
    // recomputation of the full evolution operator if dt hasn't changed.
    // Note that in d=1 we aim for performance rather than memory efficiency
    // because having kx2 and U as cached versions has O(N) memory.

    std::iota(kx2.begin(), kx2.end(), -NN / 2 + 1);
    std::rotate(kx2.begin(), kx2.begin() + NN / 2 - 1, kx2.end());
    kx2 = box.dk * box.dk * kx2 * kx2;
}

// Limit max phase change to pi/2
double Kinetic::next_dt(const SimState &state) const { return M_PI / k2_max; }

void Kinetic::step(SimState &state, const double dt) const {
    state.transform(SimState::Representation::Momentum);

    // There is no point in evolving modes with no significant power
    const double psi_k_min = 1e-12;
    for (int i = 0; i <= box.N / 2; ++i) {
        if (std::abs(state.psi[i]) < psi_k_min) {
            k2_max = kx2[i];
            break;
        }
    }

    // We only need to update the matrix exponential, if an altered time step
    // size is required
    if (dt - dt_last != 0) { // FP compare is ok since we assign dt_last to dt
        U = exp(-0.5i * kx2 * dt);
    }
    state.psi *= U;
    state.tau_aux += dt;

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

} // namespace Schroedinger
