#ifndef __SCHROEDINGER_PC_CAYLEY__
#define __SCHROEDINGER_PC_CAYLEY__

#include "driver.h"

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/SymmetricMatrix.h>

// TODO Add integrator description

namespace Schroedinger {

class PCCayley : public DefaultDriver<PCCayley> {
    const Cosmology& cosmo;
    size_t N;
    double dx;
    double dt;
    double a;
    std::unique_ptr<Interaction> potential;
    blaze::SymmetricMatrix<blaze::CompressedMatrix<double>> K;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> psi_old;
    blaze::DynamicVector<double> V_old;
    blaze::DynamicVector<std::complex<double>> dl, d, du;
    blaze::DynamicVector<std::complex<double>> du2;
    blaze::DynamicVector<int> ipiv;

   public:
    PCCayley(const Parameters& p, const SimState& state,
             const Cosmology& cosmo_);
    void step(SimState& state, const double dt);
    double next_dt(const SimState& state) const;
    REGISTER(PCCayley)
};

}  // namespace Schroedinger
#endif
