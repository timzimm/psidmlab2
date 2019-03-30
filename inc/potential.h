#ifndef __POTENTIAL_INTERFACE__
#define __POTENTIAL_INTERFACE__
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <complex>
#include <memory>
#include <utility>
#include "common.h"
#include "fftw3.h"

// Forward Declarations
struct Parameters;

// In the current version of the code performance is not our main objective,
// modularity, however, is. Hence, we implement solvers via dynamic
// polymorphism.
//
// TODO: Switch to CRTP/something template based once we decided on a solver.
namespace Potential {

class Algorithm {
   public:
    // TODO Do something useful with the second paramter or find a better
    // interface
    virtual void operator()(SimState& state) = 0;
    virtual ~Algorithm() = default;
};

namespace Poisson {

// Solves Poisson Equation by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.
// Assumes periodic boundary conditions.
class FFT : public Algorithm {
    using RCV = blaze::DynamicVector<double>;
    using RRV = blaze::DynamicVector<double, blaze::rowVector>;
    using CCV = blaze::DynamicVector<std::complex<double>>;
    using RDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<double>>;
    size_t N;
    double L;

    CCV fft;
    RDM inv_k_sq;

    fftw_plan forwards;
    fftw_plan backwards;

   public:
    constexpr static double epsilon = 1e-8;

    FFT(const Parameters& p);
    ~FFT();
    void operator()(SimState& state) override;
};
}  // namespace Poisson

// Add new potential algorithms here
enum class AlgorithmType { Poisson_FFT };

// Potential solver factory. Add new solvers here.
struct Factory {
    static std::unique_ptr<Algorithm> create(const AlgorithmType& t,
                                             const Parameters& p) {
        switch (t) {
            case AlgorithmType::Poisson_FFT:
                return std::make_unique<Poisson::FFT>(p);
            default:
                return nullptr;
        }
    }
};
}  // namespace Potential

#endif
