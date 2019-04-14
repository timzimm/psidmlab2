#ifndef __INTERFACES__
#define __INTERFACES__
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

class Observable : public Factory<Observable, Parameters> {
   public:
    virtual void compute(SimState& state) = 0;
    virtual ~Observable() = default;
};
#endif
