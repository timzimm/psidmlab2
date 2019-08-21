#ifndef __STRANG_SPLITTING__
#define __STRANG_SPLITTING__

#include "interfaces.h"

namespace Schroedinger {

template <typename FlowMapA, typename FlowMapB>
class StrangSplitting : public SchroedingerMethod::Registrar<
                            StrangSplitting<FlowMapA, FlowMapB>> {
    // Subproblem Integrators (A, B possible nonlinear)
    FlowMapA phiA;  // integrator for i del_t psi = A psi
    FlowMapB phiB;  // integrator for i del_t psi = B psi

   public:
    StrangSplitting(const Parameters& p, const SimState& state,
                    const Cosmology& cosmo_);

    // Next time step for total problem determined by stability of the
    // subproblems A and B, i.e. taking the minium of both
    double next_dt(const SimState& state) const;

    // Strang type evolution operator for i del_t psi = (A + B) psi
    void step(SimState& state, const double dt);
};

}  // namespace Schroedinger
#endif
