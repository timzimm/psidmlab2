#include <cassert>

#include "cosmology.h"
#include "parameters.h"
#include "schroedinger/uso_dkd.h"

namespace Schroedinger {

USO_DKD::USO_DKD(const Parameters& p, const SimState& state,
                 const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      k_squared(N),
      forwards(nullptr),
      backwards(nullptr) {
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.end() - (N + 1) / 2,
                k_squared.end());
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;
    auto in = const_cast<fftw_complex*>(
        reinterpret_cast<const fftw_complex*>(state.psis.data()));

    int dist_state = state.psis.spacing();

    forwards =
        fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state, in,
                           nullptr, 1, dist_state, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards = fftw_plan_many_dft(1, &N, state.M, in, nullptr, 1, dist_state,
                                   in, nullptr, 1, dist_state, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
}

USO_DKD::~USO_DKD() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}
void USO_DKD::step_internal(SimState& state, const double dt) {
    using namespace blaze;
    CCM& psis = state.psis;
    const double& a_t = state.a;
    const double a_t_dt = cosmo.a_of_tau(state.tau + dt);

    // Drift step
    psis = expand(exp(-0.5 * cmplx(0, 1) * 0.5 * (a_t + a_t_dt) * state.V * dt),
                  state.M) %
           psis;

    auto state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(forwards, state_p, state_p);

    auto kick = expand(exp(-0.5 * cmplx(0, 1) * k_squared * dt), state.M);
    psis = 1.0 / N * kick % psis;

    state_p = reinterpret_cast<fftw_complex*>(psis.data());
    fftw_execute_dft(backwards, state_p, state_p);

    pot->solve(state);

    psis = expand(exp(-0.5 * cmplx(0, 1) * 0.5 * (a_t + a_t_dt) * state.V * dt),
                  state.M) %
           psis;

    // state is now @ tau + dt
    pot->solve(state);
    state.tau += dt;
    state.a = a_t_dt;
}

void USO_DKD::step(SimState& state) {
    using namespace blaze;
    CCM& psis = state.psis;
    const double tau_final = state.tau + state.dtau;
    const double a_final = cosmo.a_of_tau(tau_final);
    const double& t = state.tau;

    auto next_dt = [&]() {
        return std::min({L * L / (2 * N * N * M_PI),
                         M_PI / (2 * a_final * max(abs(state.V)))});
    };

    for (double dt = next_dt(); t + dt < tau_final; dt = next_dt())
        step_internal(state, dt);

    step_internal(state, tau_final - t);

    state.n += 1;
}
}  // namespace Schroedinger
