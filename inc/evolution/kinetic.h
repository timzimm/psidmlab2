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
    // Wavenumbers per dimension. This is used to represent the k-grid without
    // allocating memory of O(Nx*Ny*Nz). In 2,3D the required memory is neglible
    // compared to sizeof(state). In 1D it is of the same size, tough.
    blaze::DynamicVector<double> kx;
    blaze::DynamicVector<double> ky;
    blaze::DynamicVector<double> kz;

   public:
    Kinetic(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);
};
}  // namespace Schroedinger
#endif
