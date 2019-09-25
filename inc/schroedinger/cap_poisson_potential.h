#ifndef __SCHROEDINGER_CAP_POISSON__
#define __SCHROEDINGER_CAP_POISSON__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

// Analytical time evolution operator for the free Schroedinger equation:
//          i del_t psi = -1/2 del^2_x psi
// assuming periodic boundary conditions.

class CAPPoissonPotential : public DefaultDriver<CAPPoissonPotential> {
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    const int N;                           // No. of spatial points
    const double L;                        // size of domain
    const double dx;
    const double strength;
    double dt_last;
    blaze::DynamicVector<double, blaze::columnVector> CAP;
    blaze::DynamicVector<double, blaze::columnVector> attenuator;

   public:
    CAPPoissonPotential(const Parameters& p, const SimState& state,
                        const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
