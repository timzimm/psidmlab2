#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include "fftw3.h"
#include "interfaces.h"

// Gross Pitaevskii type non-linearity of the form |psi|^2*psi

class GP : public Interaction {
   public:
    GP(const Parameters& p);
    ~GP();
    void solve(SimState& state) override;
    void solve(blaze::DynamicVector<double>& V,
               const blaze::DynamicVector<double>& source) override;
    REGISTER(GP)
};

#endif
