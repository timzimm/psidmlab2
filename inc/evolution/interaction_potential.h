#ifndef __SCHROEDINGER_INTERACTION_POTENTIAL__
#define __SCHROEDINGER_INTERACTION_POTENTIAL__

#include "domain.h"
#include "driver.h"

namespace Schroedinger {

// Exact time evolution operator for the Schroedinger equation:
//        i ∂_t psi(x,t) =   a(t')  U[|psi(x,t)|^2] psi(x,t)
// for a real parameter t'
// we set t' =  state.tau_aux, which the caller sets at its own will.

class InteractionPotential : public DefaultDriver<InteractionPotential> {
    const Cosmology& cosmo;
    const Domain box;
    std::unique_ptr<Interaction> pot;

   public:
    InteractionPotential(const Parameters& p, const SimState& state,
                         const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(InteractionPotential)
};
}  // namespace Schroedinger
#endif
