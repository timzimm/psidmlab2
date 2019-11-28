#ifndef __POISSON_FD__
#define __POISSON_FD__

#include "fftw3.h"
#include "interfaces.h"

namespace Poisson {

// Solves Poisson Equation in 1D by approximating the laplacian via a second
// order finite difference approximation. The resulting cyclic tridiagonal
// matrix is then inverted by relying on a LU-decomposition implemented by
// LAPACK.
class FD : public Interaction {
    using RCV = blaze::DynamicVector<double>;
    using RCM = blaze::DynamicMatrix<double, blaze::columnMajor>;
    double dx;
    size_t N;
    // LU-decompostion of poisson matrix which never changes
    RCV dl, d, du, du2;
    blaze::DynamicVector<int> ipiv;  // pivot array

   public:
    FD(const Parameters& p);
    void solve(SimState& state) override;
    void solve(blaze::DynamicVector<double, blaze::columnVector>& V,
               const blaze::DynamicVector<double, blaze::columnVector>& source)
        override;
    REGISTER(FD)
};

}  // namespace Poisson
#endif
