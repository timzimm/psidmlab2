#ifndef __SCHROEDINGER_PC_CAYLEY__
#define __SCHROEDINGER_PC_CAYLEY__

#include "driver.h"

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/SymmetricMatrix.h>

// TODO Add integrator description

namespace Schroedinger {

class PCCayley : public DefaultDriver<PCCayley> {
    // integer row vector
    using IRV = blaze::DynamicVector<int>;
    // complex row vector
    using CRV = blaze::DynamicVector<std::complex<double>>;
    // real column vector
    using RCV = blaze::DynamicVector<double, blaze::columnVector>;
    // dense complex matrix (mixed state)
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;
    // sparse symmetric matrix (second derivative approcimation)
    using RSM = blaze::SymmetricMatrix<blaze::CompressedMatrix<double>>;

    const Cosmology& cosmo;
    size_t N;
    size_t M;
    double dx;
    double dt;
    double a;
    std::unique_ptr<Interaction> potential;
    RSM K;
    CCM psi_old;
    RCV V_old;
    CRV dl, d, du;
    CRV du2;
    IRV ipiv;

   public:
    PCCayley(const Parameters& p, const SimState& state,
             const Cosmology& cosmo_);
    void step(SimState& state, const double dt);
    double next_dt(const SimState& state) const;
    REGISTER(PCCayley)
};

}  // namespace Schroedinger
#endif
