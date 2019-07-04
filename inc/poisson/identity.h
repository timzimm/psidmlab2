#ifndef __POISSON_IDENTITY__
#define __POISSON_IDENTITY__

#include "interfaces.h"

// Identity for the static potential case.
// Essentially only needed to avoid touching the code in main

namespace Poisson {

class Identity : public PotentialMethod::Registrar<Identity> {
   public:
    Identity(const Parameters& p);
    void solve(SimState& state) override;
    void solve(blaze::DynamicVector<double, blaze::columnVector>& V,
               const blaze::DynamicVector<double, blaze::columnVector>& source)
        override;
};

}  // namespace Poisson
#endif
