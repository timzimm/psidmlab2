#ifndef __SCREENED_POISSON_FC__
#define __SCREENED_POISSON_FC__

#include "blaze/math/dense/DynamicVector.h"
#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Approximately computes
//
//  V[psi](x) = = ∫ dx' -1/(2*eps) * exp(-eps*|x-x'|) |psi(x')|^2
//  with x in [box.xmin, box.xmax]
//
// by assuming |psi|^2 decays fast (exponentially) to zero so that
//
// (i) Dirichlet conditions are applicable for psi (and only psi)
// (ii) the integration domain can be truncated to [box.xmin,box.xmax].
//
// Expressing |psi|^2(x) via a sin-pseudospactral representation
// yields integrals over the domain which can be solved analytically.

namespace ScreenedPoisson {

class FastConvolution : public Interaction {
    const Domain box;
    const double epsilon;
    blaze::DynamicVector<double> source_k;
    fftw_plan_ptr dst;

   public:
    FastConvolution(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(FastConvolution)
};

}  // namespace ScreenedPoisson
#endif
