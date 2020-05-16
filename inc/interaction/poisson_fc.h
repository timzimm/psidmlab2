#ifndef __POISSON_FC__
#define __POISSON_FC__

#include "blaze/math/dense/DynamicVector.h"
#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Approximately computes
//
//  V[psi](x) = = ∫ dx' 1/2 * |x-x'| |psi(x')|^2    x in [box.xmin, box.xmax]
//
// by assuming |psi|^2 decays fast (exponentially) to zero so that
//
// (i) Dirichlet conditions are applicable for psi (and only psi)
// (ii) the integration domain can be truncated to [box.xmin,box.xmax].
//
// Expressing |psi|^2(x) via a sin-pseudospactral representation
// yields integrals over the domain which can be solved analytically.

namespace Poisson {

class FastConvolution : public Interaction {
    const Domain box;
    blaze::DynamicVector<double> source_k;
    fftw_plan_ptr dst;

   public:
    FastConvolution(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(FastConvolution)
};

}  // namespace Poisson
#endif
