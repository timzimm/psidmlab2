#include "evolution/pc_cayley.h"
#include "blaze_utils.h"
#include "cosmology.h"
#include "parameters.h"

#include <blaze/math/Column.h>

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

PCCayley::PCCayley(const Parameters& p, const SimState& state,
                   const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      N{p["Simulation"]["N"].get<size_t>()},
      dx{p["Simulation"]["L"].get<double>() / N},
      dt{p["Simulation"]["dtau"].get<double>()},
      a{cosmo.a_of_tau(-1)},
      potential{Interaction::make(
          p["Simulation"]["potential"].get<std::string>(), p)},
      K(N, N),
      psi_old(N, 1),
      V_old{state.V},
      dl(N),
      d(N),
      du(N),
      du2(N),
      ipiv(N) {
    blaze::CompressedMatrix<double> dx2_coeff(N, N);
    dx2_coeff.reserve(3 * N);
    dx2_coeff.append(0, 0, -2);
    dx2_coeff.append(0, 1, 1);
    dx2_coeff.append(0, N - 1, 1);
    dx2_coeff.finalize(0);

    for (int i = 1; i < N - 1; ++i) {
        dx2_coeff.append(i, i - 1, 1);
        dx2_coeff.append(i, i, -2);
        dx2_coeff.append(i, i + 1, 1);
        dx2_coeff.finalize(i);
    }
    dx2_coeff.append(N - 1, 0, 1);
    dx2_coeff.append(N - 1, N - 2, 1);
    dx2_coeff.append(N - 1, N - 1, -2);
    dx2_coeff.finalize(N - 1);
    K = declsym(-1.0 / (2 * dx * dx) * dx2_coeff);
}

void PCCayley::step(SimState& state, const double dt) {
    const std::complex<double> I{0, 1.0};
    const double a_next = cosmo.a_of_tau(state.tau + dt);

    state.transform(SimState::Representation::Position);
    potential->solve(state);
    // Save the input state for the corrector step - no copy
    swap(column(psi_old, 0), state.psi);
    swap(V_old, state.V);

    const auto a_V = a * V_old;

    // PREDICTOR STEP
    state.psi = psi_old - 0.5i * dt * (K * psi_old + a_V * psi_old);

    // Left hand side matrix M+
    dl = du = -1.0i / (4 * dx * dx) * dt;
    d = -2.0 * dl + 0.5i * dt * a * V_old + 1.0;

    // Solve cyclic tridiagonal matrix equation
    gctrf(dl, d, du, du2, ipiv);
    gctrs(dl, d, du, du2, ipiv, expand(state.psi, 1));

    // CORRECTOR STEP
    potential->solve(state);
    auto V_corr = 0.5 * (a_V + a_next * state.V);

    state.psi = psi_old - 0.5i * dt * (K * psi_old + V_corr * psi_old);

    // Left hand side matrix M+
    // LU decomposition works in place. Thus, everything has to be
    // reinitialized.
    dl = du = -1.0i / (4 * dx * dx) * dt;
    d = -2.0 * dl + 0.5i * dt * 0.5 * (a * V_old + a_next * state.V) + 1.0,

    gctrf(dl, d, du, du2, ipiv);
    gctrs(dl, d, du, du2, ipiv, state.psi);

    state.tau += dt;
    a = a_next;
    state.n += 1;
}

// TODO Stability Analysis
double PCCayley::next_dt(const SimState& state) const { return dt; }
}  // namespace Schroedinger
