#ifndef __SCHROEDINGER_USO_DKD__
#define __SCHROEDINGER_USO_DKD__

#include <fftw3.h>

#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
#include "driver.h"

// Unitary split operator. Propagates psi according to
// psi_(n+1) = exp(-i/2 * V) exp(-i * K) exp(-i/2 * V) psi_(n)
// for separable Hamiltionian H = K + V = -1/2*del^2_x + a(t)*V(psi)

namespace Schroedinger {

class USO_DKD : public DefaultDriver<USO_DKD> {
    // Real column vector
    using RCV = blaze::DynamicVector<double, blaze::columnVector>;
    // Complex dense matrix in column-major order
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;

    // Data members
    const Cosmology& cosmo;
    std::unique_ptr<Interaction> pot;
    const int N;
    const double L;
    double dt_last;
    double a;
    RCV k_squared;
    CCM kick;

   public:
    USO_DKD(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    void step(SimState& state, const double dt);
    double next_dt(const SimState& state) const;
    REGISTER(USO_DKD);
};
}  // namespace Schroedinger
#endif
