#ifndef __STRANG_SPLITTING__
#define __STRANG_SPLITTING__

#include "interfaces.h"

template <typename FlowMapA, typename FlowMapB>
class StrangSplitting
    : public Stepper::Registrar<StrangSplitting<FlowMapA, FlowMapB>> {
    // Subproblem Integrators (A, B possible nonlinear)
    FlowMapA phiA;  // integrator for i del_t psi = A psi
    FlowMapB phiB;  // integrator for i del_t psi = B psi

   public:
    StrangSplitting(const Parameters& p, const SimState& state,
                    const Cosmology& cosmo_)
        : phiA(p, state, cosmo_), phiB(p, state, cosmo_) {}
    // Next time step for total problem determined by stability of the
    // subproblems A and B, i.e. taking the minium of both
    double next_dt(const SimState& state) const {
        return std::min({phiA.next_dt(state), phiB.next_dt(state)});
    }

    // Strang type evolution operator for i del_t psi = (A + B) psi
    void step(SimState& state, const double dt) {
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
};

#endif
