#ifndef __STATE__
#define __STATE__
#include "fftw.h"
#include "parameters_fwd.h"

#include <blaze/math/DynamicVector.h>
#include <complex>
#include <cstddef>

struct SimState {
    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    unsigned int n;  // time step number
    double tau;      // current time
    blaze::DynamicVector<double> V;

    blaze::DynamicVector<std::complex<double>> psi;

    SimState();
    void transform(const Representation target);

   private:
    Representation representation;
    fftw_plan_ptr position_to_momentum;
    fftw_plan_ptr momentum_to_position;
    int N_plan;
    double *psi_ptr;
};

void operator>>(const SimState &state, Parameters &p);

inline decltype(auto) delta_from(const SimState &state) {
    return real(conj(state.psi) * state.psi) - 1;
}

#endif
