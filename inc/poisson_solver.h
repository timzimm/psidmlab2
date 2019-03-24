#ifndef __POISSON_SOLVER_INTERFACE__
#define __POISSON_SOLVER_INTERFACE__
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "common.h"
#include "fftw3.h"

// In the current version of the code performace is not our main objective,
// modularity, however, is. Hence, we implement solvers via dynamic
// polymorphism. The ABC below defines the common interface of all solvers.
//
// TODO: Switch to CRTP once we decided on a solver

namespace Solvers {
class Poisson {
   public:
    virtual blaze::DynamicVector<double, blaze::rowVector> solve(
        blaze::DynamicVector<double, blaze::rowVector>& source) = 0;
    virtual ~Poisson() = default;
};

class FFT : public Poisson {
    using RRV = blaze::DynamicVector<double, blaze::rowVector>;
    using CCV = blaze::DynamicVector<std::complex<double>>;
    using RDM = blaze::DiagonalMatrix<blaze::DynamicMatrix<double>>;
    size_t N;
    double L;
    RRV potential;

    CCV fft;
    RDM inv_k_sq;

    fftw_plan forwards;
    fftw_plan backwards;

   public:
    constexpr static double epsilon = 1e-8;

    FFT(const Parameters& p);
    ~FFT();
    RRV solve(RRV& source) override;
};
}  // namespace Solvers
#endif
