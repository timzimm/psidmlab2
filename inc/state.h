#ifndef __STATE__
#define __STATE__
#include "parameters_fwd.h"

#include <blaze/math/DynamicVector.h>
#include <fftw3.h>
#include <complex>

struct SimState {
    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    int n;       // time step number
    double tau;  // current time
    blaze::DynamicVector<double> V;

    blaze::DynamicVector<std::complex<double>> psi;

    SimState(const Parameters& p);
    void transform(const Representation target);
    ~SimState();

   private:
    Representation representation;
    fftw_plan position_to_momentum;
    fftw_plan momentum_to_position;
    int N_plan;
};

void operator>>(const SimState& state, Parameters& p);

inline decltype(auto) delta_from(const SimState& state) {
    return real(conj(state.psi) * state.psi) - 1;
}

#endif
