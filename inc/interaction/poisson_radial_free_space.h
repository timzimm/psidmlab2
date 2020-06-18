#ifndef __POISSON_RADIAL_FREE_SPACE__
#define __POISSON_RADIAL_FREE_SPACE__

#include "blaze/math/dense/DynamicVector.h"
#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Approximately solves the FREE SPACE PROBLEM:
//                  r * ∂^2_r W = source(r)     r in [0, ∞)
// with
//      W(0) = 0
//      lim_∞ W(r) = 1
// // by TRUNCATING onto [0, box.L] and homogenizing the potential
//
//      U(r) = W(r) - r/box.L = 4πr V(r) - r/box.L
//
// Assuming HOMOGENOUS DIRICHLET CONDITIONS on the boundaries we have:
//
//                  r * ∂^2_r U = source(r)     r in [0, R]
// with
//      U(0) = U(R) = 0
//
// The solve methods return
//      V(r) = 1/(4πr) * U(r) + 1/(4π box.L)

namespace Poisson {

class RadialFreeSpace : public Interaction {
    const Domain box;
    blaze::DynamicVector<double> source_k;
    fftw_plan_ptr dst;

   public:
    RadialFreeSpace(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(RadialFreeSpace)
};

}  // namespace Poisson
#endif
