#ifndef __INTERFACES__
#define __INTERFACES__
#include "factory.h"
#include "state.h"

#include <boost/variant.hpp>

// This header defines all interfaces of
//
// Stepper (like Schroedinger::USO_DKD, Schroedinger::PCCaley, ...)
// Driver (like SimpleDriver, StableDriver, ...)
// PotentialMethod (like Poisson::FFT, Poisson::FD, ...)
// ObservableFunctor (like Observable::DensityConstrast, Observable::Potential)

// Forward Declarations
class Cosmology;
#include "parameters_fwd.h"

// Class representing a potential generator (e.g. a fixed potential,
// GP-nonlinearity) or a solver (e.g. for Poisson's equation)
class PotentialMethod : public Factory<PotentialMethod, const Parameters&> {
   public:
    virtual void solve(SimState& state) = 0;
    // Verbose interface
    virtual void solve(
        blaze::DynamicVector<double, blaze::columnVector>& V,
        const blaze::DynamicVector<double, blaze::columnVector>& source) = 0;
    virtual ~PotentialMethod() = default;
};

// Class representing a PDE stepper for Schroedingers equation (e.g. linear,
// general non-linear ...)
class Stepper : public Factory<Stepper, const Parameters&, const SimState&,
                               const Cosmology&> {
   public:
    // Step to t + dt without considering any stability criterion. It is the
    // responsibility of the derived Stepper to change the
    // representation of state.psis as required, that is no guarantees are made
    // whether the passed in state is in momentum, position or any other
    // representation.
    virtual void step(SimState& state, const double dt) = 0;
    // Returns next stable dt based on current state
    virtual double next_dt(const SimState& state) const = 0;
    virtual ~Stepper() = default;
};

class Driver : public Factory<Driver, const Parameters&> {
   public:
    // Integrates the numerical solution over the finite, macroscopic time
    // interval [state.tau, state.tau + t_final]
    virtual void integrate(std::unique_ptr<Stepper>& stepper, SimState& state,
                           const double t_final) = 0;
    virtual ~Driver() = default;
};

// Class representing an Observable constructed from the simulation state.
class ObservableFunctor
    : public Factory<ObservableFunctor, const Parameters&, const Cosmology&> {
   public:
    using ReturnType = boost::variant<
        const blaze::DynamicMatrix<double>&,
        const blaze::DynamicMatrix<double, blaze::columnMajor>&,
        const blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>&,
        const blaze::DynamicVector<double>&,
        const blaze::DynamicVector<double, blaze::columnVector>&>;

    virtual ReturnType compute(const SimState& state) = 0;
    virtual ~ObservableFunctor() = default;
};
#endif
