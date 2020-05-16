#ifndef __NONLOCAL_GP__
#define __NONLOCAL_GP__

#include "domain.h"
#include "fftw.h"
#include "interaction/periodic_convolution.h"
#include "interfaces.h"

// Computes the interaction:
//
// V[psi](x) = G * |psi|^2 + gamma*|psi|^2  x in [box.xmin, box.xmax]
//
// with periodic boundary conditions for psi and V.
//
// Convolution is solved by delegating to PeriodicConvolution

class NonlocalGP : public Interaction {
    PeriodicConvolution poisson;
    double gamma;

   public:
    NonlocalGP(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(NonlocalGP)
};

#endif
