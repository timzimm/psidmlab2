#ifndef __POISSON_FFT_EPSILON__
#define __POISSON_FFT_EPSILON__

#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Solves Poisson Equation by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.
// Assumes periodic boundary conditions.

namespace Poisson {

class FFTEpsilon : public Interaction {
    const Domain box;
    const double epsilon;
    fftw_plan_ptr fwd;
    fftw_plan_ptr bwd;
    double *real_ptr;

   public:
    FFTEpsilon(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(FFTEpsilon)
};

}  // namespace Poisson
#endif
