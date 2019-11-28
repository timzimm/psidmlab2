#ifndef __INTERFACES__
#define __INTERFACES__
#include "factory.h"
#include "state.h"

#include <boost/variant.hpp>

// This header defines all interfaces of
//
// Interaction (like Schroedinger::USO_DKD, Schroedinger::PCCaley, ...)
// Driver (like SimpleDriver, StableDriver, ...)
// PotentialMethod (like Poisson::FFT, Poisson::FD, ...)
// ObservableFunctor (like Observable::DensityConstrast, Observable::Potential)

// Forward Declarations
class Cosmology;
#include "parameters_fwd.h"

// Class representing the generic operator U[psi] which could be ...
// U[psi] = V_extern(x)                         (linear SE)
// U[psi] = g|psi|^2(x)                         (GPE)
// U[psi] = V(x) with d^2_x V(x) = |psi|^2 - 1  (SPE)
// U[psi] = G * (|psi|^2 - 1)                   (nonlocal NLSE)
class Interaction : public Factory<Interaction, const Parameters&> {
   public:
    using InterfaceType = Interaction;
    // given the wavefunction state.psi compute state.V
    virtual void solve(SimState& state) = 0;
    // Verbose interface with same functionality.
    // Output is written into vector V
    virtual void solve(
        blaze::DynamicVector<double, blaze::columnVector>& V,
        const blaze::DynamicVector<double, blaze::columnVector>& source) = 0;
    virtual ~Interaction() = default;
};

// Class representing a time evolution operator of the Schroedinger equation
class TimeEvolution : public Factory<TimeEvolution, const Parameters&,
                                     const SimState&, const Cosmology&> {
   public:
    using InterfaceType = TimeEvolution;
    // Step to t + dt without considering any stability criterion. It is the
    // responsibility of the derived Stepper to change the
    // representation of state.psis as required, that is no guarantees are made
    // whether the passed in state is in momentum, position or any other
    // representation.
    virtual void step(SimState& state, const double dt) = 0;
    // Returns next stable dt based on current state
    virtual double next_dt(const SimState& state) const = 0;
    // Integrates the numerical solution over the finite, macroscopic time
    // interval [state.tau, state.tau + t_final]
    virtual void integrate(SimState& state, const double t_final) = 0;
    virtual ~TimeEvolution() = default;
};

// Class representing an Observable constructed from the simulation state.
class ObservableFunctor
    : public Factory<ObservableFunctor, const Parameters&, const Cosmology&> {
   public:
    using InterfaceType = ObservableFunctor;
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
