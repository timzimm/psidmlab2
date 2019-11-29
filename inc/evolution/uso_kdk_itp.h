#ifndef __SCHROEDINGER_USO_KDK__
#define __SCHROEDINGER_USO_KDK__

#include "driver.h"
#include "interfaces.h"

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>

namespace Schroedinger {

class USO_KDK_ITP : public DefaultDriver<USO_KDK_ITP> {
    // Real column vector

    const Cosmology& cosmo;
    std::unique_ptr<Interaction> pot;
    blaze::DynamicVector<double> pot_external;
    const int N;
    const double L;
    const double norm0;
    blaze::DynamicVector<double> k_squared;
    blaze::DynamicVector<std::complex<double>> kick;
    double dt_last;

   public:
    USO_KDK_ITP(const Parameters& p, const SimState& state,
                const Cosmology& cosmo_);
    double next_dt(const SimState& state) const;
    void step(SimState& state, const double dt);

    REGISTER(USO_KDK_ITP)
};
}  // namespace Schroedinger

#endif
