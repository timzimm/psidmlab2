#ifndef __PERIODIC_CONVOLUTION__
#define __PERIODIC_CONVOLUTION__

#include "blaze/math/dense/DynamicVector.h"
#include "domain.h"
#include "fftw.h"
#include "interfaces.h"

// Computes the solution to the non-local interaction:
//
//        xmax                            ∞
// V(x) = ∫ G(x,x') |psi(x')|^2 dx' = 1/L ∑ G_k (|psi|^2)_k e^{2 pi i/L(x-xmin)}
//        xmin                           k=1
//
// with PERIODIC BOUNDARY CONDITIONS for psi and V:
// The series is truncated after box.N terms and computed via a DFT accelerated
// by a FFT. For this to work we assume |psi|^2 is sufficiently smooth so that
// (|psi|^2)_k decays rapidly and dominates the (asymptotic) 1/k^2 decay of all
// G_k below.
//
// The following interactions are supported:
// 1D Poisson:
//   G(x,x') = 1/2 |x-x'| - 1/2 [(x-x')^2/L + L/6]
//   G_k     = -1.0/k^2
// 1D Screened Poisson:
//   G(x,x') = -1/(2*eps) cosh(eps*(L/2 - |x-x'|))/sinh(eps*L/2) + 1/(L*eps^2)
//   G_k     = -1.0/(k^2 + eps^2)
// Periodic Line Adiabatic Model:
//   G(x,x') = ---
//   G_k = -1/(4pi) U(1, 1, 1/2*(k*eps*2)) (confluent hypergeometric function)
//
// The computation only relies on G_k!
//
// To avoid recomputation of the kernel every time the interaction is computed,
// we store it in the vector G_k
//
// Memory Requiremnets: O(N)
// Time Complexity: O(NlogN)

class PeriodicConvolution : public Interaction {
    const Domain box;
    blaze::DynamicVector<double> G_k;
    fftw_plan_ptr fwd;
    fftw_plan_ptr bwd;
    double *real_ptr;
    enum class InteractionType { Poisson, ScreenedPoisson, LAM };

   public:
    PeriodicConvolution(const Parameters &p, const SimState &state);
    void solve(SimState &state) override;
    void solve(blaze::DynamicVector<double> &V,
               const blaze::DynamicVector<double> &source) override;
    REGISTER(PeriodicConvolution)
};

#endif
