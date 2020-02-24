#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include "fftw.h"
#include "interfaces.h"
#include "state.h"

// Solves Poisson Equation by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.
// Assumes periodic boundary conditions.

namespace Poisson {

class FFT : public Interaction {
    const size_t N;
    const double L;
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
