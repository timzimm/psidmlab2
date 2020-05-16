#ifndef __DNGF_BACKWARD_FORWARD_EULER_SINE_PSEUDOSPECTRAL__
#define __DNGF_BACKWARD_FORWARD_EULER_SINE_PSEUDOSPECTRAL__

#include "domain.h"
#include "driver.h"
#include "fftw.h"

class Cosmology;

namespace DNGF {

// Approximates the Continuous Normalized Gradient Flow (CNGF) by means of a
// Forward-Backward-Euler time discrtization, in which linear terms are
// approximated forward in time and non-linear terms backward in time:
//
// (i)      1/∆t * (phi* - phi^n) =
//               1/2 ∂_xx phi* - alpha phi* + alpha*phi^n - a V^n phi^n
//
// (ii)     phi^(n+1) = phi* / || phi*||
//
// and phi(box.xmin) = phi(box.xmax) = 0 (HOMOGENOUES DIRICHLET)
// the alpha term is used as a stabilizer to increase the range of stable ∆t.
// Discretization of space is done via the sine-pseudospectral approach.
//
// See W.Bao Journal of Compuatational Physics 219 (2006) 836-854 for more
// information

class BFSP : public DefaultDriver<BFSP> {
    const Domain box;
    const Cosmology& cosmo;
    mutable blaze::DynamicVector<double> phi;
    mutable double bmin;
    mutable double bmax;
    mutable double alpha;
    const double norm0;
    fftw_plan_ptr dst;
    std::unique_ptr<Interaction> pot;

   public:
    BFSP(const Parameters& p, const SimState& state, const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(BFSP)
};
}  // namespace DNGF
#endif
