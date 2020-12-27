#ifndef __STATE__
#define __STATE__
#include "domain.h"
#include "fftw.h"

#include <blaze/math/DynamicVector.h>
#include <complex>
#include <cstddef>

struct SimState {
    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    unsigned int n; // time step number
    double tau;     // current time
    double dtau;    // last time increment
    double tau_aux; // auxillary time
    blaze::DynamicVector<double> V;

    blaze::DynamicVector<std::complex<double>> psi;

    void transform(const Representation target);
    SimState(const Domain &box);

  private:
    Representation representation;
    fftw_plan_ptr position_to_momentum;
    fftw_plan_ptr momentum_to_position;
    int N_plan;
    size_t N_transform;
    double norm;
    double *psi_ptr;
};

inline decltype(auto) delta_from(const SimState &state) {
    return real(conj(state.psi) * state.psi) - 1;
}
inline decltype(auto) rho_from(const SimState &state) {
    return real(conj(state.psi) * state.psi);
}

#endif
