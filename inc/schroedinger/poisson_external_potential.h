#ifndef __SCHROEDINGER_POISSON_EXTERNAL_POTENTIAL__
#define __SCHROEDINGER_POISSON_EXTERNAL_POTENTIAL__

#include "driver.h"

#include <fstream>

namespace Schroedinger {

class PoissonExternalPotential
    : public DefaultDriver<PoissonExternalPotential> {
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    blaze::DynamicVector<double> pot_external;

   public:
    PoissonExternalPotential(const Parameters& p, const SimState& state,
                             const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
