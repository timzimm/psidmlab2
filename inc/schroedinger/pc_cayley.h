#ifndef __SCHROEDINGER_PC_CAYLEY__
#define __SCHROEDINGER_PC_CAYLEY__

#include <complex>
#include <memory>
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "blaze/math/SymmetricMatrix.h"
#include "cosmology.h"
#include "interfaces.h"

// TODO Add integrator description

namespace Schroedinger {

class PCCayley : public SchroedingerMethod::Registrar<PCCayley> {
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
    std::unique_ptr<PotentialMethod> potential;
    RSM K;  // Cyclic Kinetic matrix (i.e. the - 0.5*second derivative)
    CCM psi_old;
    RCV V_old;
    CRV dl, d, du;  // diagonal and off diagonal elements of matrix M+
    CRV du2;   // second super diagonal of upper tridiagonal matrix in LU decomp
    IRV ipiv;  // pivot vector for lapack routines

   public:
    PCCayley(const Parameters& p, const SimState& state,
             const Cosmology& cosmo_);
    void step(SimState& state) override;
};

}  // namespace Schroedinger
#endif
