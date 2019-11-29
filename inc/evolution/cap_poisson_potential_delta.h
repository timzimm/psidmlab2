#ifndef __SCHROEDINGER_CAP_POISSON_DELTA__
#define __SCHROEDINGER_CAP_POISSON_DELTA__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

class CAPPoissonPotentialDelta
    : public DefaultDriver<CAPPoissonPotentialDelta> {
    const Cosmology& cosmo;
    std::unique_ptr<Interaction> pot;
    const int N;
    const double L;
    const double dx;
    const double strength;
    double dt_last;
    blaze::DynamicVector<double> CAP;
    blaze::DynamicVector<double> attenuator;

   public:
    CAPPoissonPotentialDelta(const Parameters& p, const SimState& state,
                             const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
