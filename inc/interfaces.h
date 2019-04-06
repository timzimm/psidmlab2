#ifndef __INTERFACES__
#define __INTERFACES__
#include "factory.h"

// This header defines all interfaces of
//
// SchroedingerMethod (like Schroedinger::USO_DKD, Schroedinger::CN, ...)
// PotentialMethod (like Poisson::FFT, Wave::LAX, ...)
//
// The interface is currently as generic as possible and might change in the
// future to something more expressive/functional.

// Forward Declarations
struct SimState;
struct Parameters;

// Note that a PotentialMethod can be based on a Stepper (wave equation)
// or a solver (poisson equation).
class PotentialMethod : public Factory<PotentialMethod, Parameters> {
   public:
    virtual void operator()(SimState& state) = 0;
    virtual ~PotentialMethod() = default;
};

class SchroedingerMethod : public Factory<SchroedingerMethod, Parameters> {
   public:
    virtual void operator()(SimState& state) = 0;
    virtual ~SchroedingerMethod() = default;
};
#endif
