#ifndef __SCHROEDINGER_INTERACTION_POTENTIAL_TRAPEZODIAL__
#define __SCHROEDINGER_INTERACTION_POTENTIAL_TRAPEZODIAL__

#include "driver.h"

#include <fstream>

namespace Schroedinger {

// (A) Exact time evolution operator for the Schroedinger equation:
//        i del_t psi(x,t) =  [ a U[psi(x,t)] + V_ext(x) ] psi(x,t)
// for a = const. and
// (B) 2nd order approximate evolution operator for non-autonomous
// Hamiltonian:
//        i del_t psi(x,t) =  [ a(t) U[psi(x,t)] + V_ext(x) ] psi(x,t)

class InteractionPotentialTrapezodial
    : public DefaultDriver<InteractionPotentialTrapezodial> {
    const Cosmology& cosmo;
    double dt_last;
    int N;
    std::unique_ptr<Interaction> pot;
    blaze::DynamicVector<double> pot_external;

   public:
    InteractionPotentialTrapezodial(const Parameters& p, const SimState& state,
                                    const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
