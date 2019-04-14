#ifndef __INTERFACES__
#define __INTERFACES__
#include "blaze/math/DynamicMatrix.h"
#include "factory.h"

// This header defines all interfaces of
//
// SchroedingerMethod (like Schroedinger::USO_DKD, Schroedinger::CN, ...)
// PotentialMethod (like Poisson::FFT, Poisson::FD, ...)
//
// The interface is currently as generic as possible and might change in the
// future to something more expressive/functional.

// Forward Declarations
#define __FORWARD__
struct SimState;
#include "parameters_fwd.h"

class PotentialMethod : public Factory<PotentialMethod, Parameters> {
   public:
    virtual void solve(SimState& state) = 0;
    virtual ~PotentialMethod() = default;
};

class SchroedingerMethod : public Factory<SchroedingerMethod, Parameters> {
   public:
    virtual void step(SimState& state) = 0;
    virtual ~SchroedingerMethod() = default;
};

class ObservableFunctor : public Factory<ObservableFunctor, Parameters> {
   public:
    // virtual function templates do not exists. So double is the way to go.
    virtual blaze::DynamicMatrix<double, blaze::columnMajor> compute(
        const SimState& state) = 0;
    virtual ~ObservableFunctor() = default;
};
#endif
