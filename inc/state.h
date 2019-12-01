#ifndef __STATE__
#define __STATE__
#include "parameters_fwd.h"

#include <blaze/math/DynamicVector.h>
#include <fftw3.h>
#include <array>
#include <complex>

/*
 * SimState holds the discretized version of all fields relevant for the
 * integration. Currently this includes the pure state wavefunction psi(x) as
 * well as the gravitaional potential V(x). The type of the discretization is
 * controlled by the SimState::Representation enum and currently comprises of:
 *
 * SimState::Representation::Position: psi on a uniform grid with dx = L/dof
 * SimState::Representation::Momentum: psi in momentum basis with dk = 2pi/L
 *
 * The field discrizations of SimState are populated by the ICGenerator class.
 * The rest of its memebers are essentially the parameters of the "Domain
 * Properties" section in the config file.
 */

struct SimState {
    using interval_t = std::array<double, 2ul>;

    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    size_t n;                     // time step number
    double tau;                   // current time
    std::vector<int> dofs;        // degrees of freedom per dimension
    std::vector<interval_t> box;  // box limits
    size_t dims;                  // # of dimensions
    size_t N_total;               // total # of grid points (N_1 * ... * N_dims)

    // Fields
    blaze::DynamicVector<double> V;
    blaze::DynamicVector<std::complex<double>> psi;

    SimState(const Parameters& p);
    void transform(const Representation target);
    ~SimState();

   private:
    Representation representation;
    fftw_plan position_to_momentum;
    fftw_plan momentum_to_position;
};

void operator>>(const SimState& state, Parameters& p);

inline decltype(auto) delta_from(const SimState& state) {
    return real(conj(state.psi) * state.psi) - 1;
}

#endif
