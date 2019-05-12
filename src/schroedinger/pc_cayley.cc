#include "schroedinger/pc_cayley.h"
#include "blaze/math/UniformVector.h"
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
    // Compute K at @ construction (symmetry automatically enforced)
    blaze::band(K, -1) = blaze::UniformVector<double>(N - 1, 1);
    blaze::band(K, 0) = blaze::UniformVector<double>(N, -2);
    K(N - 1, 0) = 1.0;
    K *= -1.0 / (2 * dx * dx);
}

void PCCayley::step(SimState& state) {
    using namespace blaze;

    std::complex<double> I{0, 1.0};
    double dt = state.dtau;
    double t = state.tau;
    double a = state.a;
    double a_da = cosmo.a_of_tau(t + dt);

    // Save the input state for the corrector step - no copy
    swap(psi_old, state.psis);
    swap(V_old, state.V);

    auto a_V = a * expand(V_old, M);

    // PREDICTOR STEP
    state.psis = psi_old - I * 0.5 * dt * (K * psi_old + a_V % psi_old);

    // Left hand side matrix M+
    dl = du = -1.0 * std::complex<double>{0, 1} / (4 * dx * dx) * dt;
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
    // TODO: Can we change this even though LAPACK defines this behavior?
    dl = du = -1.0 * std::complex<double>{0, 1} / (4 * dx * dx) * dt;
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
