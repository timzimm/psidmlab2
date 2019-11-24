#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include "fftw3.h"
#include "interfaces.h"

// Gross Pitaevskii type non-linearity of the form |psi|^2*psi

namespace Poisson {

class GP : public PotentialMethod {
   public:
    GP(const Parameters& p);
    ~GP();
    void solve(SimState& state) override;
    void solve(blaze::DynamicVector<double, blaze::columnVector>& V,
               const blaze::DynamicVector<double, blaze::columnVector>& source)
        override;
    REGISTER(GP)
};

}  // namespace Poisson
#endif
