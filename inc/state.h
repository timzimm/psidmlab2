#ifndef __STATE__
#define __STATE__
#include "domain.h"
#include "fftw.h"
#include "parameters_fwd.h"

#include <blaze/math/DynamicVector.h>
#include <complex>

// Forward Declaration

/*
 * SimState holds the discretized version of all fields relevant for the
 * integration. Currently this includes the pure state wavefunction psi(x) as
 * well as the gravitational potential V(x). The type of the discretization is
 * controlled by the SimState::Representation enum and currently comprises of:
 *
 * SimState::Representation::Position: psi on a uniform grid with dx = L/dof
 * SimState::Representation::Momentum: psi in momentum basis with dk = 2pi/L
 *
 */

struct SimState {
    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    const Domain& box;
    size_t n;    // time step number
    double tau;  // current time

    // Fields
    blaze::DynamicVector<double> V;
    blaze::DynamicVector<std::complex<double>> psi;

    SimState(const Domain& p);
    void transform(const Representation target);

   private:
    Representation representation;
    fftw_plan_ptr position_to_momentum;
    fftw_plan_ptr momentum_to_position;
};

inline decltype(auto) delta_from(const SimState& state) {
    return real(conj(state.psi) * state.psi) - 1;
}

#endif
