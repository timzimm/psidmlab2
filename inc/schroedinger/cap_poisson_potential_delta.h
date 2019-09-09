#ifndef __SCHROEDINGER_CAP_POISSON_DELTA__
#define __SCHROEDINGER_CAP_POISSON_DELTA__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

// Analytical time evolution operator for the free Schroedinger equation:
//          i del_t psi = -1/2 del^2_x psi
// assuming periodic boundary conditions.

class CAPPoissonPotentialDelta
    : public DefaultDriver<CAPPoissonPotentialDelta> {
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    const int N;                           // No. of spatial points
    const double L;                        // size of domain
    const double dx;
    const double strength;
    const double width;
    double dt_last;
    blaze::DynamicVector<double, blaze::columnVector> CAP;
    blaze::DynamicVector<double, blaze::columnVector> attenuator;

   public:
    CAPPoissonPotentialDelta(const Parameters& p, const SimState& state,
                             const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
