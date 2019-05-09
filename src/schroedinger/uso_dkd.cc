#include "schroedinger/uso_dkd.h"
#include <cassert>
#include "cosmology.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

USO_DKD::USO_DKD(const Parameters& p, const SimState& state,
                 const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      K(N, N),
      D(N, N),
      k_squared(N),
      forwards(nullptr),
      backwards(nullptr) {
    for (int k = 0; k < N / 2; ++k) {
        k_squared[k] = 2 * M_PI / L * k;
        k_squared[k] *= k_squared[k];
    }
    for (int k = N / 2; k < N; ++k) {
        k_squared[k] = 2 * M_PI / L * (k - N);
        k_squared[k] *= k_squared[k];
    }
    // This makes me sad. On so many levels :(
    const auto const_in =
        reinterpret_cast<const fftw_complex*>(state.psis.data());
    auto in = const_cast<fftw_complex*>(const_in);

    int dist_state = state.psis.spacing();

    forwards =
        fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state, in,
                           nullptr, 1, dist_state, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards = fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state,
                                   in, nullptr, 1, dist_state, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
}

// Kick Operator - returns a blaze expression
decltype(auto) USO_DKD::kick(const CCM& psis_in_k, const double dt,
                             const double weight) {
    auto diag_K = blaze::diagonal(K);
    diag_K = blaze::exp(-1.0 * weight / 2 * cmplx(0, 1) * k_squared * dt);

    return K * psis_in_k;
}

// Drift Operator - returns a blaze expression
decltype(auto) USO_DKD::drift(const CCM& psis_in_x, const RCV& V,
                              const double t, const double dt,
                              const double weight) {
    auto diag_D = blaze::diagonal(D);
    diag_D =
        blaze::exp(-1.0 * weight * cmplx(0, 1) * cosmo.a_of_tau(t) * V * dt);
    // This is a dense-sparse product
    return D * psis_in_x;
}

USO_DKD::~USO_DKD() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}

void USO_DKD::step(SimState& state) {
    double dt = state.dtau;
    double t = state.tau;
    CCM& psis = state.psis;
    RCV& V = state.V;

    // Drift step
    // It is crucial to always use the most up-to date version of psi.
    // Therefore, we cannot reuse the potential of the last step but
    // have to recalculate it.
    pot->solve(state);
    psis = blaze::evaluate(drift(psis, V, t, dt, 1.0 / 2));

    fftw_execute(forwards);

    psis = blaze::evaluate(1.0 / N * kick(psis, dt, 1.0));

    fftw_execute(backwards);

    pot->solve(state);
    psis = blaze::evaluate(drift(psis, V, t, dt, 1.0 / 2));

    // psis and V are now @ tau + dtau
    // Update time information
    state.tau += dt;
    state.n += 1;
    state.a = cosmo.a_of_tau(t + dt);
}
}  // namespace Schroedinger
