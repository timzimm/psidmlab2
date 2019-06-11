#include "schroedinger/uso_kdk.h"
#include <cassert>
#include "cosmology.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

USO_KDK::USO_KDK(const Parameters& p, const SimState& state,
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
      backwards(nullptr),
      forwards_op(nullptr),
      backwards_op(nullptr),
      psis_cached(N, state.M) {
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    for (int k = 0; k < N / 2; ++k) {
        k_squared[k] = 2 * M_PI / L * k;
        k_squared[k] *= k_squared[k];
    }
    for (int k = N / 2; k < N; ++k) {
        k_squared[k] = 2 * M_PI / L * (k - N);
        k_squared[k] *= k_squared[k];
    }

    const auto const_in =
        reinterpret_cast<const fftw_complex*>(state.psis.data());
    auto in = const_cast<fftw_complex*>(const_in);
    auto out = reinterpret_cast<fftw_complex*>(psis_cached.data());

    int dist_state = state.psis.spacing();
    int dist_cache = psis_cached.spacing();

    // in-place
    forwards =
        fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state, in,
                           nullptr, 1, dist_state, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards = fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state,
                                   in, nullptr, 1, dist_state, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);

    // out-of-place
    forwards_op = fftw_plan_many_dft(
        1, &N, state.M, in, nullptr, 1, state.psis.spacing(), out, nullptr, 1,
        psis_cached.spacing(), FFTW_FORWARD, FFTW_ESTIMATE);
    backwards_op = fftw_plan_many_dft(
        1, &N, state.M, out, nullptr, 1, psis_cached.spacing(), in, nullptr, 1,
        state.psis.spacing(), FFTW_BACKWARD, FFTW_ESTIMATE);

    // Init cached wavefunction in k
    fftw_execute(forwards_op);

    // Normalize
    psis_cached /= N;
}

// Kick Operator
decltype(auto) USO_KDK::kick(const CCM& psis_in_k, const double dt,
                             const double weight) {
    auto diag_K = blaze::diagonal(K);

    diag_K = blaze::exp(-1.0 * weight / 2 * cmplx(0, 1) * k_squared * dt);

    // This is a dense-sparse product
    return K * psis_in_k;
}

// Drift Operator
decltype(auto) USO_KDK::drift(const CCM& psis_in_x, const RCV& V,
                              const double t, const double dt,
                              const double weight) {
    auto diag_D = blaze::diagonal(D);

    // TODO Fix a(tau)
    diag_D =
        blaze::exp(-1.0 * weight * cmplx(0, 1) * cosmo.a_of_tau(t) * V * dt);

    // This is a dense-sparse product
    return D * psis_in_x;
}

USO_KDK::~USO_KDK() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
    fftw_destroy_plan(forwards_op);
    fftw_destroy_plan(backwards_op);
}

void USO_KDK::step(SimState& state) {
    const double dt = state.dtau;
    const double t = state.tau;
    CCM& psis = state.psis;
    RCV& V = state.V;

    // We can spare the initial FFT if we use the cached and normalized psi_in_k
    // representation of the last step
    psis = blaze::evaluate(kick(psis_cached, dt, 1.0 / 2));

    fftw_execute(backwards);

    // Recompute potential
    pot->solve(state);
    psis = blaze::evaluate(drift(psis, V, t, dt, 1.0));

    fftw_execute(forwards);

    // Update cached psis and normalize
    psis_cached = blaze::evaluate(1.0 / N * kick(psis, dt, 1.0 / 2));

    fftw_execute(backwards_op);

    // psis and V are now @ tau + dtau
    state.tau += dt;
    state.n += 1;
    state.a = cosmo.a_of_tau(t + dt);
}

}  // namespace Schroedinger
