#ifndef __SCHROEDINGER_INTERACTION_EXTERNAL_POTENTIAL__
#define __SCHROEDINGER_INTERACTION_EXTERNAL_POTENTIAL__

#include "driver.h"

#include <fstream>

namespace Schroedinger {

class InteractionExternalPotential
    : public DefaultDriver<InteractionExternalPotential> {
    const Cosmology& cosmo;  // Cosmological model for a(t)
    double dt_last;
    int N;
    std::unique_ptr<PotentialMethod> pot;  // potential method
    blaze::DynamicVector<double> pot_external;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>
        drift;  // evolution operator in x space

   public:
    InteractionExternalPotential(const Parameters& p, const SimState& state,
                                 const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
