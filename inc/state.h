#ifndef __STATE__
#define __STATE__
#include "parameters_fwd.h"

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <fftw3.h>
#include <complex>

struct SimState {
    // Representation tags for the stored wavefunction
    enum class Representation { Position, Momentum };

    int n;       // time step number
    double tau;  // current time
    blaze::DynamicVector<double> V;

    // state = sum_i lambda_i * |psi_i><psi_i|
    int M;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> psis;
    blaze::DynamicVector<double> lambda;

    SimState(const Parameters& p);
    void transform(const Representation target);
    ~SimState();

   private:
    Representation representation;
    fftw_plan position_to_momentum;
    fftw_plan momentum_to_position;
    int N_plan, M_plan;
};

void operator>>(const SimState& state, Parameters& p);

inline decltype(auto) delta_from(const SimState& state) {
    auto psi2 = real(conj(state.psis) % state.psis) - 1.0;
    return psi2 * state.lambda;
}

#endif
