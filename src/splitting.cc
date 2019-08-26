.include "splitting.h"

    template <typename FlowMapA, typename FlowMapB>
    StrangSplitting<FlowMapA, FlowMapB>::StrangSplitting(
        const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : phiA(p, state, cosmo_),
    phiB(p, state, cosmo_) {}

template <typename FlowMapA, typename FlowMapB>
double StrangSplitting<FlowMapA, FlowMapB>::next_dt(
    const SimState& state) const {
    return std::min({phiA.next_dt(state), phiB.next_dt(state)});
}

template <typename FlowMapA, typename FlowMapB>
void StrangSplitting<FlowMapA, FlowMapB>::step(SimState& state,
                                               const double dt) {
    const double t = state.tau;
    const int n = state.n;

    phiA.step(state, 0.5 * dt);
    state.tau = t;
    state.n = n;
    phiB.step(state, dt);
    state.tau = t;
    state.n = n;
    phiA.step(state, 0.5 * dt);
    state.tau = t;
    state.n = n;

    state.tau = t + dt;
    state.n += 1;
}
