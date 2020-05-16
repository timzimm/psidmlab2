#ifndef __DNGF_BACKWARD_FORWARD_EULER_FOURIER_PSEUDOSPECTRAL__
#define __DNGF_BACKWARD_FORWARD_EULER_FOURIER_PSEUDOSPECTRAL__

#include <memory>
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
// and PERIODIC BOUNDARY CONDITIONS on phi.
// the alpha term is used as a stabilizer to increase the range of stable ∆t.
// Discretization of space is done via the fourier-pseudospectral approach.
//
// See W.Bao Journal of Compuatational Physics 219 (2006) 836-854 for more
// information

class BFFP : public DefaultDriver<BFFP> {
    const Domain box;
    const Cosmology& cosmo;
    mutable blaze::DynamicVector<std::complex<double>> aVphi_k;
    mutable double bmin;
    mutable double bmax;
    mutable double alpha;
    const double norm0;
    std::unique_ptr<Interaction> pot;
    fftw_plan_ptr r2c;
    fftw_plan_ptr c2r;

   public:
    BFFP(const Parameters& p, const SimState& state, const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt) const;

    REGISTER(BFFP)
};
}  // namespace DNGF
#endif
