#ifndef __STATE__
#define __STATE__
#include "parameters_fwd.h"

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <complex>

struct SimState {
    int n;        // time step number
    double tau;   // current time
    double dtau;  // current time increment
    double a;     // current scale factor
    blaze::DynamicVector<double> V;

    // state = sum_i lambda_i * |psi_i><psi_i|
    int M;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor> psis;
    blaze::DynamicVector<double> lambda;

    SimState(const Parameters& p);
};

void operator>>(const SimState& state, Parameters& p);

inline decltype(auto) delta_from(const SimState& state) {
    auto psi2 = real(conj(state.psis) % state.psis) - 1;
    return psi2 * state.lambda;
}

#endif
