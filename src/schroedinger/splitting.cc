#include "schroedinger/splitting.h"

namespace Schroedinger {

template <typename FlowMapA, typename FlowMapB>
StrangSplitting<FlowMapA, FlowMapB>::StrangSplitting(const Parameters& p,
                                                     const SimState& state,
                                                     const Cosmology& cosmo_)
    : phiA(p, state, cosmo_), phiB(p, state, cosmo_) {}

template <typename FlowMapA, typename FlowMapB>
double StrangSplitting<FlowMapA, FlowMapB>::next_dt(
    const SimState& state) const {
    return std::min({phiA.next_dt(), phiB.next_dt()});
}

template <typename FlowMapA, typename FlowMapB>
void StrangSplitting<FlowMapA, FlowMapB>::step(SimState& state,
                                               const double dt) {
    const double t = state.tau;
    const int n = state.n;

    phiA.step(state, 0.5 * dt);
    phiB.step(state, dt);
    phiA.step(state, 0.5 * dt);

    state.tau = t + dt;
    state.n += 1;
}

/* template class StrangSplitting<Kinetic, PeriodicPotential>; */
/* template class StrangSplitting<Kinetic, PMLPeriodicPotential>; */

}  // namespace Schroedinger

