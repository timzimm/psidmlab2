#include "schroedinger/pc_cayley.h"
#include "blaze_utils.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

PCCayley::PCCayley(const Parameters& p, const SimState& state,
                   const Cosmology& cosmo_)
    : cosmo{cosmo_},
      N{p["Simulation"]["N"].get<size_t>()},
      M{p["Initial Conditions"]["M"].get<size_t>()},
      dx{p["Simulation"]["L"].get<double>() / N},
      potential{PotentialMethod::make(
          p["Simulation"]["potential"].get<std::string>(), p)},
      K(N, N),
      psi_old(N, M),
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

void PCCayley::step(SimState& state) {
    using namespace blaze;

    const std::complex<double> I{0, 1.0};
    const double dt = state.dtau;
    const double t = state.tau;
    const double a = state.a;
    const double a_da = cosmo.a_of_tau(t + dt);

    // Save the input state for the corrector step - no copy
    swap(psi_old, state.psis);
    swap(V_old, state.V);

    const auto a_V = a * expand(V_old, M);

    // PREDICTOR STEP
    state.psis = psi_old - I * 0.5 * dt * (K * psi_old + a_V % psi_old);

    // Left hand side matrix M+
    dl = du = -1.0 * I / (4 * dx * dx) * dt;
    d = -2.0 * dl + I * 0.5 * dt * a * V_old + 1.0;

    // Solve cyclic tridiagonal matrix equation
    gctrf(dl, d, du, du2, ipiv);
    gctrs(dl, d, du, du2, ipiv, state.psis);

    // CORRECTOR STEP
    potential->solve(state);
    auto V_corr = 0.5 * (a_V + a_da * expand(state.V, M));

    state.psis = psi_old - I * 0.5 * dt * (K * psi_old + V_corr % psi_old);

    // Left hand side matrix M+
    // LU decomposition works in place. Thus, everything has to be
    // reinitialized.
    dl = du = -1.0 * I / (4 * dx * dx) * dt;
    d = -2.0 * dl + I * 0.5 * dt * 0.5 * (a * V_old + a_da * state.V) + 1.0,

    gctrf(dl, d, du, du2, ipiv);
    gctrs(dl, d, du, du2, ipiv, state.psis);

    // At last we calculate the potential again such that state is @ t + dt
    potential->solve(state);
    state.tau = t + dt;
    state.n += 1;
    state.a = a_da;
}

}  // namespace Schroedinger
