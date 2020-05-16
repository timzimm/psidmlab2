#ifndef __SCHROEDINGER_INTERACTION_POTENTIAL_MAGNUS__
#define __SCHROEDINGER_INTERACTION_POTENTIAL_MAGNUS__

#include "domain.h"
#include "driver.h"

namespace Schroedinger {

// (A) Exact time evolution operator for the Schroedinger equation:
//        i ∂_t psi(x,t) =   a U[|psi(x,t)|^2] psi(x,t)
// for a = const. and
//
// (B) 2nd order approximate evolution operator for non-autonomous
// Hamiltonian:
//        i ∂_t psi(x,t) =  a(t) U[|psi(x,t)|^2] psi(x,t)
// by approximating the time evolution operator up to first order in Magnus:
//
// U(t,t0) = exp( -i ∫dt a(t) U[|psi(x,t)|^2] ) + O(∆t^3)
//         = exp( -i U[|psi(x,t0)|^2] ∫dt a(t) ) + O(∆t^3)
//         ≈ exp( -i U[|psi(t0)|^2] a(t + ∆t/2) ∆t ) + O(∆t^3)
//
// Note that |psi(t)|^2 is a conserved quantity for both (A) and (B). Thus,
// taking U out of the integral is NOT an approximation but in fact EXACT.

class InteractionPotentialMagnus
    : public DefaultDriver<InteractionPotentialMagnus> {
    const Cosmology& cosmo;
    const Domain box;
    std::unique_ptr<Interaction> pot;

   public:
    InteractionPotentialMagnus(const Parameters& p, const SimState& state,
                               const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(InteractionPotentialMagnus)
};
}  // namespace Schroedinger
#endif
