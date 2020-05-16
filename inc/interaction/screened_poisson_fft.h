#ifndef __YUKAWA_FFT__
#define __YUKAWA_FFT__

#include "domain.h"
#include "fftw.h"
#include "interfaces.h"
// Computes the solution to the screened Poisson equation:
//
//       ∂_xx V -  epsilon^2 V =  |psi|^2 - 1   x in [box.xmin, box.xmax]
//
// with periodic boundary conditions for psi and V.
//
// The Equation is solved by taking a discrete Fourier transformation of
// the source, transforming its coefficient and taking its inverse FT.

namespace ScreenedPoisson {

class FFT : public Interaction {
    const Domain box;
    const double epsilon;
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

}  // namespace ScreenedPoisson
#endif
