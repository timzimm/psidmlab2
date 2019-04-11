#ifndef __POISSON_FD__
#define __POISSON_FD__

#include <complex>
#include "blaze/math/DynamicVector.h"
#include "fftw3.h"
#include "interfaces.h"

// Solves Poisson Equation by approximatating the laplacian via a second order
// finite difference approximation. The resulting cyclic tridiagonal matrix
// is then inverted by relying on a LU-decomposition implemented by LAPACK.

namespace Poisson {

class FD : public PotentialMethod::Registrar<FD> {
    using RCV = blaze::DynamicVector<double>;
    double dx;
    size_t N;
    // LU-decompostion of poisson matrix which never changes
    RCV dl, d, du, du2;
    blaze::DynamicVector<int> ipiv;  // pivot array

   public:
    FD(const Parameters& p);
    ~FD();
    void solve(SimState& state) override;
};

}  // namespace Poisson
#endif
