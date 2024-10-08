#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Computes the solution to:
//
//                 ∂_xx V =  (|psi|^2 - 1)    x in [box.xmin, box.xmax]
//
// with periodic boundary conditions for psi and V.
//
// Poisson Equation is solved by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.

namespace Poisson {

class FFT : public Interaction {
    const Domain box;
    fftw_plan_ptr fwd;
    fftw_plan_ptr bwd;
    double *real_ptr;

   public:
    FFT(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(FFT)
};

}  // namespace Poisson
#endif
