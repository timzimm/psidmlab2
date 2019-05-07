#ifndef __INTERFACES__
#define __INTERFACES__
#include <boost/variant.hpp>
#include <complex>
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "factory.h"

// This header defines all interfaces of
//
// SchroedingerMethod (like Schroedinger::USO_DKD, Schroedinger::CN, ...)
// PotentialMethod (like Poisson::FFT, Poisson::FD, ...)
// ObservableFunctor (like Observable::DensityConstrast, Observable::Potential)
//
// The interface is currently as generic as possible and might change in the
// future to something more expressive/functional.

// Forward Declarations
#define __FORWARD__
struct SimState;
class Cosmology;
#include "parameters_fwd.h"

class PotentialMethod : public Factory<PotentialMethod, const Parameters&> {
   public:
    virtual void solve(SimState& state) = 0;
    virtual ~PotentialMethod() = default;
};

class SchroedingerMethod : public Factory<SchroedingerMethod, const Parameters&,
                                          const SimState&, const Cosmology&> {
   public:
    virtual void step(SimState& state) = 0;
    virtual ~SchroedingerMethod() = default;
};

class ObservableFunctor : public Factory<ObservableFunctor, const Parameters&> {
   public:
    using ReturnType = boost::variant<
        const blaze::DynamicMatrix<double, blaze::columnMajor>&,
        const blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>&,
        const blaze::DynamicVector<double>&>;

    virtual ReturnType compute(const SimState& state) = 0;
    virtual ~ObservableFunctor() = default;
};
#endif
