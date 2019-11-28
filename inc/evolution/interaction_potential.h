#ifndef __SCHROEDINGER_INTERACTION_POTENTIAL__
#define __SCHROEDINGER_INTERACTION_POTENTIAL__

#include "driver.h"

namespace Schroedinger {

class InteractionPotential : public DefaultDriver<InteractionPotential> {
    const Cosmology& cosmo;
    std::unique_ptr<Interaction> pot;

   public:
    InteractionPotential(const Parameters& p, const SimState& state,
                         const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
