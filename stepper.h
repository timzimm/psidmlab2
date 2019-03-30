#ifndef __POISSON_SOLVER_INTERFACE__
#define __POISSON_SOLVER_INTERFACE__
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <complex>
#include <memory>
#include <utility>
#include "fftw3.h"

// Forward Declarations
struct Parameters;

// In the current version of the code performance is not our main objective,
// modularity, however, is. Hence, we implement solvers via dynamic
// polymorphism. The ABC below defines the common interface of all solvers,
// independent of the underlying elliptic PDE to solve.
//
// TODO: Switch to CRTP/something template based once we decided on a solver.

template <typename T>
class Stepper {
   public:
    virtual std::pair<blaze::DynamicVector<T, blaze::rowVector>, double>
    operator()(const blaze::DynamicVector<T, blaze::rowVector>& u,
               const double dt) = 0;
    virtual ~Solver() = default;
};

#include "stepper/uso_dkd.h"

namespace Schroedinger {

// Add new solver types here
enum class Type { FFT };

// Solver factory. Add new solvers here.
struct Factory {
    static std::unique_ptr<::Solver> create(const Type& t,
                                            const Parameters& p) {
        switch (t) {
            case Type::FFT:
                return std::make_unique<FFT>(p);
            default:
                return nullptr;
        }
    }
};

}  // namespace Schroedinger

#endif
