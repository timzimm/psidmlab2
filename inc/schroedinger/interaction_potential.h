#ifndef __SCHROEDINGER_INTERACTION_POTENTIAL__
#define __SCHROEDINGER_INTERACTION_POTENTIAL__

#include "driver.h"

namespace Schroedinger {

// Analytical time evolution operator for the free Schroedinger equation:
//          i del_t psi = -1/2 del^2_x psi
// assuming periodic boundary conditions.

class InteractionPotential : public DefaultDriver<InteractionPotential> {
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method

   public:
    InteractionPotential(const Parameters& p, const SimState& state,
                         const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
