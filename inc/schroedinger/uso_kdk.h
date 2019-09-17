#ifndef __SCHROEDINGER_USO_KDK__
#define __SCHROEDINGER_USO_KDK__

#include "driver.h"
#include "interfaces.h"

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>

// Unitary split operator. Propagates psi according to
// psi_(n+1) = exp(-i/2 * K) exp(-i * V) exp(-i/2 * K) psi_(n)
// for separable Hamiltionian H = K + V = -1/2*del^2_x + a(t)*V(psi)

namespace Schroedinger {

class USO_KDK : public DefaultDriver<USO_KDK> {
    // Real column vector
    using RCV = blaze::DynamicVector<double, blaze::columnVector>;
    // Complex dense matrix in column-major order
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;

    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    const int N;                           // No. of spatial points
    const double L;                        // size of domain
    RCV k_squared;
    CCM kick;
    double dt_last;

   public:
    USO_KDK(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);

    REGISTER(USO_KDK)
};
}  // namespace Schroedinger

#endif
