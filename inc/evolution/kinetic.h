#ifndef __SCHROEDINGER_KINETIC__
#define __SCHROEDINGER_KINETIC__

#include "domain.h"
#include "driver.h"

class Cosmology;

namespace Schroedinger {

// Exact time evolution operator for the equation:
//
//                  i ∂_t psi = -1/2 ∂^2_x psi
//                    ∂_t v   = 1
//
// assuming PERIODIC BOUNDARY CONDITIONS. The time evolution operator reads:
//
// U(t, ∆t) (psi) = ( ƒ^{-1} exp(-i/2 * k^2 * ∆t) ƒ )
//          (v)     (          v + ∆t               )
//
// with ƒ as DCT operator.

class Kinetic : public DefaultDriver<Kinetic> {
    const Domain box;
    mutable double dt_last;
    mutable double k2_max;
    blaze::DynamicVector<double> kx2;
    mutable blaze::DynamicVector<std::complex<double>> U;

  public:
    Kinetic(const Parameters &p, const SimState &state,
            const Cosmology &cosmo_);
    double next_dt(const SimState &state) const;
    void step(SimState &state, const double dt) const;

    REGISTER(Kinetic)
};
} // namespace Schroedinger
#endif
