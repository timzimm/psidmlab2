#include "schroedinger/uso_kdk.h"
#include "cosmology.h"
#include "parameters.h"

#include <algorithm>

namespace Schroedinger {

USO_KDK::USO_KDK(const Parameters& p, const SimState& state,
                 const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
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
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;

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

USO_KDK::~USO_KDK() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
    fftw_destroy_plan(forwards_op);
    fftw_destroy_plan(backwards_op);
}

void USO_KDK::step_internal(SimState& state, const double dt) {
    using namespace blaze;
    CCM& psis = state.psis;
    const double& a_t = state.a;
    const double a_t_dt = cosmo.a_of_tau(state.tau + dt);

    // Set kick operator once for the entire time step
    auto kick = expand(exp(-0.5 * 0.5 * cmplx(0, 1) * k_squared * dt), state.M);

    // We can spare the initial FFT if we use the cached and normalized
    // psi_in_k representation of the last step
    psis = kick % psis_cached;

    auto state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(backwards, state_p, state_p);

    // Recompute potential
    pot->solve(state);

    // Set drift operator with intermediate potential
    auto drift = expand(
        exp(-1.0 * cmplx(0, 1) * 0.5 * (a_t + a_t_dt) * state.V * dt), state.M);
    psis = drift % psis;

    state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(forwards, state_p, state_p);

    // Update cached psis and normalize
    psis_cached = 1.0 / N * kick % psis;

    // psis is now @ tau + dt
    // potential is in intermediate state
    state.tau += dt;
    state.a = a_t_dt;
}

void USO_KDK::step(SimState& state) {
    using namespace blaze;
    CCM& psis = state.psis;
    const double tau_final = state.tau + state.dtau;
    const double a_final = cosmo.a_of_tau(tau_final);
    const double& t = state.tau;

    auto next_dt = [&]() {
        return std::min(
            {L * L / (N * N * M_PI), M_PI / (a_final * max(abs(state.V)))});
    };

    for (double dt = next_dt(); t + dt < tau_final; dt = next_dt())
        step_internal(state, dt);

    step_internal(state, tau_final - t);

    auto state_p = reinterpret_cast<fftw_complex*>(psis.data());
    auto cache_p = reinterpret_cast<fftw_complex*>(psis_cached.data());
    fftw_execute_dft(backwards_op, cache_p, state_p);

    state.n += 1;
}

}  // namespace Schroedinger
