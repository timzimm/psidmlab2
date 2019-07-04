#include "schroedinger/uso_kdk.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>

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
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.end() - (N + 1) / 2,
                k_squared.end());
    k_squared *= 2 * M_PI / L;
    k_squared *= k_squared;

    auto in = const_cast<fftw_complex*>(
        reinterpret_cast<const fftw_complex*>(state.psis.data()));
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
    forwards_op =
        fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state, out,
                           nullptr, 1, dist_cache, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards_op = fftw_plan_many_dft(1, &N, state.M, out, nullptr, 1,
                                      dist_cache, in, nullptr, 1, dist_state,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);

    // Init cached wavefunction in k
    fftw_execute(forwards_op);

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
    using namespace blaze;
    const double dt = state.dtau;
    const double t = state.tau;
    CCM& psis = state.psis;
    RCV& V = state.V;

    // Set kick operator once for the entire time step
    auto diag_K = blaze::diagonal(K);
    diag_K = blaze::exp(-0.5 * 0.5 * cmplx(0, 1) * k_squared * dt);

    // We can spare the initial FFT if we use the cached and normalized psi_in_k
    // representation of the last step
    psis = K * psis_cached;

    auto state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(backwards, state_p, state_p);

    // Recompute potential
    pot->solve(state);

    // Set drift operator with intermediate potential
    auto diag_D = blaze::diagonal(D);
    diag_D = blaze::exp(-1.0 * cmplx(0, 1) * cosmo.a_of_tau(t) * V * dt);

    psis = D * psis;

    state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(forwards, state_p, state_p);

    // Update cached psis and normalize
    psis_cached = 1.0 / N * K * psis;

    auto cache_p = reinterpret_cast<fftw_complex*>(psis_cached.data());
    fftw_execute_dft(backwards_op, cache_p, state_p);

    // psis and V are now @ tau + dtau
    state.tau += dt;
    state.n += 1;
    state.a = cosmo.a_of_tau(t + dt);
}

}  // namespace Schroedinger
