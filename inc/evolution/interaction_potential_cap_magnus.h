#ifndef __SCHROEDINGER_CAP_MAGNUS__
#define __SCHROEDINGER_CAP_MAGNUS__

#include "domain.h"
#include "driver.h"

class Cosmology;

namespace Schroedinger {

// 2nd order approximate evolution operator for non-autonomous
// Hamiltonian:
// i ∂_t psi(x,t) =  a(t) U[|psi(x,t)|^2] psi(x,t) - i/2 phi(x) psi(x,t)

class InteractionPotentialCAPMagnus
    : public DefaultDriver<InteractionPotentialCAPMagnus> {
    const Cosmology& cosmo;
    const Domain box;
    mutable double dt_last;
    std::unique_ptr<Interaction> pot;
    blaze::DynamicVector<double> phi;
    mutable blaze::DynamicVector<double> attenuator;

   public:
    InteractionPotentialCAPMagnus(const Parameters& p, const SimState& state,
                                  const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(InteractionPotentialCAPMagnus);
};
}  // namespace Schroedinger
#endif
