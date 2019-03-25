#ifndef __POISSON_SOLVER_INTERFACE__
#define __POISSON_SOLVER_INTERFACE__
#include <memory>
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "common.h"
#include "fftw3.h"

namespace Poisson {

// In the current version of the code performance is not our main objective,
// modularity, however, is. Hence, we implement solvers via dynamic
// polymorphism. The ABC below defines the common interface of all solvers.
//
// TODO: Switch to CRTP once we decided on a solver.
class Solver {
   public:
    virtual blaze::DynamicVector<double, blaze::rowVector> solve(
        const blaze::DynamicVector<double, blaze::rowVector>& source) = 0;
    virtual ~Solver() = default;
};

class FFT : public Solver {
    using RRV = blaze::DynamicVector<double, blaze::rowVector>;
    using CCV = blaze::DynamicVector<std::complex<double>>;
    using RDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<double>>;
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
    RRV solve(const RRV& source) override;
};

// Add new solver types here
enum class Type { FFT };

// Solver factory. Add new solvers here.
struct Factory {
    static std::unique_ptr<Solver> create(const Type& t, const Parameters& p) {
        switch (t) {
            case Type::FFT:
                return std::make_unique<FFT>(p);
            default:
                return nullptr;
        }
    }
};
}  // namespace Poisson
#endif
