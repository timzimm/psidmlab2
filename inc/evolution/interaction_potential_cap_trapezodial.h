#ifndef __SCHROEDINGER_CAP_INTERACTION__
#define __SCHROEDINGER_CAP_INTERACTION__

#include "domain.h"
#include "driver.h"

class Cosmology;

namespace Schroedinger {

// 2nd order approximate evolution operator for non-autonomous
// Hamiltonian:
// i ∂_t psi(x,t) =  a(t) U[|psi(x,t)|^2] psi(x,t) - i/2 phi(x) psi(x,t)

class InteractionPotentialCAPTrapezodial
    : public DefaultDriver<InteractionPotentialCAPTrapezodial> {
    const Cosmology& cosmo;
    const Domain box;
    std::unique_ptr<Interaction> pot;
    mutable double dt_last;
    blaze::DynamicVector<double> CAP;
    mutable blaze::DynamicVector<double> attenuator;

   public:
    InteractionPotentialCAPTrapezodial(const Parameters& p,
                                       const SimState& state,
                                       const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(InteractionPotentialCAPTrapezodial);
};
}  // namespace Schroedinger
#endif
