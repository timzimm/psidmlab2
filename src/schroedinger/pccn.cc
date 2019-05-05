#include "schroedinger/pccn.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/UniformVector.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

PCCN::PCCN(const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : cosmo{cosmo_},
      dx{p["Simulation"]["L"].get<double>() / p["Simulation"]["N"].get<int>()},
      N{p["Simulation"]["N"].get<size_t>()},
      potential{PotentialMethod::make(
          p["Simulation"]["potential"].get<std::string>(), p)},
      K(N, N) {
    // Cyclic Kinetic matrix never changes - compute it @ construction
    blaze::band(K, -1) = blaze::UniformVector<double>(N - 1, 1);
    blaze::band(K, 0) = blaze::UniformVector<double>(N, -2);
    blaze::band(K, 1) = blaze::UniformVector<double>(N - 1, 1);
    K(0, N - 1) = 1.0;
    K(N - 1, 0) = 1.0;

    K *= 1.0 / (2 * dx * dx);
}

void PCCN::step(SimState& state) {
    using namespace blaze;
    // sparse real diagonal matrix
    using RDM = DiagonalMatrix<CompressedMatrix<double>>;

    std::complex<double> I{0, 1.0};
    double dt = state.dtau;
    double t = state.tau;
    double a = state.a;
    double a_da = cosmo.a_of_tau(t + dt);

    // Make explicit copies, since we need the original state in the corrector
    // step;
    CCM psis = state.psis;
    RDM V(N, N);
    diagonal(V) = state.V;

    // Predictor step - Forward Euler
    state.psis = psis + dt * I * (K * psis - a * V * psis);

    // Corrector step
    potential->solve(state);
    RDM V_pred(N, N);
    diagonal(V_pred) = state.V;

    state.psis = psis + 0.5 * dt * I * K * (psis + state.psis);
    state.psis -= 0.5 * I * (a * V * psis + a_da * V_pred * state.psis);

    // At last we calculate the potential again such that state is @ t + dt
    potential->solve(state);
    state.tau = t + dt;
    state.n += 1;
    state.a = a_da;
}

}  // namespace Schroedinger
