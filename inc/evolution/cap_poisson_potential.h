#ifndef __SCHROEDINGER_CAP_POISSON__
#define __SCHROEDINGER_CAP_POISSON__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

class CAPPoissonPotential : public DefaultDriver<CAPPoissonPotential> {
    const Cosmology& cosmo;
    std::unique_ptr<Interaction> pot;
    const int N;
    const double L;
    const double dx;
    const double strength;
    const double width;
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
