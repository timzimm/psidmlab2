#ifndef __SCHROEDINGER_KINETIC__
#define __SCHROEDINGER_KINETIC__

#include "driver.h"

class Cosmology;

namespace Schroedinger {

// Exact time evolution operator for the free Schroedinger equation:
//          i del_t psi = -1/2 del^2_x psi
// assuming periodic boundary conditions.

class Kinetic : public DefaultDriver<Kinetic> {
    double dt_last;
    const int N;
    const double L;
    blaze::DynamicVector<double> kx2;

   public:
    Kinetic(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
