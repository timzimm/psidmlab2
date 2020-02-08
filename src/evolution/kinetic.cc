#include "evolution/kinetic.h"
#include "parameters.h"

#include <algorithm>

namespace Schroedinger {
using namespace blaze;
using namespace std::complex_literals;

Kinetic::Kinetic(const Parameters& p, const SimState& state, const Cosmology&)
    : DefaultDriver(p),
      dt_last{-1},
      kx(state.box.Ns[0], 0),
      ky(state.box.Ns[1], 0),
      kz(state.box.Ns[2], 0) {
    const int Nx = state.box.Ns[0];
    const int Ny = state.box.Ns[1];
    const int Nz = state.box.Ns[2];

    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    //          f0 f1 ... fN/2, f-N/2+1 ... f-1
    // TODO: This ordering us unfortunate. Without it we could use linspace()
    // and get away without any memory allocation. Maybe a more general vector
    // generator can do the trick.
    std::iota(kx.begin(), kx.end(), -1 * Nx / 2 + 1);
    std::rotate(kx.begin(), kx.begin() + Nx / 2 - 1, kx.end());

    // Special treatment is required for Y and Z because these axes might
    // be dummy dimensions (Ni=0). In k-space only k=0 should exist in this
    // case.
    auto next_even = [](int i) {
        return static_cast<int>(std::round(i / 2.0) * 2.0);
    };
    std::iota(ky.begin(), ky.end(), -1 * next_even(Ny) / 2 + 1);
    std::rotate(ky.begin(), ky.begin() + next_even(Ny) / 2 - 1, ky.end());
    std::iota(kz.begin(), kz.end(), -1 * next_even(Nz) / 2 + 1);
    std::rotate(kz.begin(), kz.begin() + next_even(Nz) / 2 - 1, kz.end());

    // Frequency increment per dimension.
    kx *= 2 * M_PI / state.box.box_lengths[0];
    ky *= 2 * M_PI / state.box.box_lengths[1];
    kz *= 2 * M_PI / state.box.box_lengths[2];
}

// Limit max phase change to pi/2 per step
double Kinetic::next_dt(const SimState& state) const {
    return std::max({std::pow(state.box.box_lengths[0], 2) /
                         (std::pow(state.box.Ns[0], 2) * M_PI),
                     std::pow(state.box.box_lengths[1], 2) /
                         (std::pow(state.box.Ns[1], 2) * M_PI),
                     std::pow(state.box.box_lengths[2], 2) /
                         (std::pow(state.box.Ns[2], 2) * M_PI)});
}

void Kinetic::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Momentum);

    auto ex = uniform(state.box.Ns[0], 1.0);
    auto ey = uniform(state.box.Ns[1], 1.0);
    auto ez = uniform(state.box.Ns[2], 1.0);

    // The entire k^2-grid can than be represented as sum of Kronecker
    // products.
    auto k_squared = kron(kron(ez, ey), kx * kx) + kron(kron(ez, ky * ky), ex) +
                     kron(kron(kz * kz, ey), ex);

    // In principle we could reuse exp(...) if dt=dt_last.
    // However, this requires again storage of O(Nx*Ny*Nz) which is
    // unacceptable. Hence we compute the exponential every time.
    state.psi *= exp(-0.5i * k_squared * dt);

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
