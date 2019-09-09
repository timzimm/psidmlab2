#ifndef __SCHROEDINGER_KINETIC__
#define __SCHROEDINGER_KINETIC__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

// Analytical time evolution operator for the free Schroedinger equation:
//          i del_t psi = -1/2 del^2_x psi
// assuming periodic boundary conditions.

class Kinetic : public DefaultDriver<Kinetic> {
    double dt_last;  // holds dt of last step call
    const int N;     // No. of spatial points
    const double L;  // size of domain
    blaze::DynamicVector<double, blaze::columnVector> k_squared;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>
        kick;  // evolution operator in k space

   public:
    Kinetic(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
