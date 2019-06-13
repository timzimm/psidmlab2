#ifndef __STATE__
#define __STATE__
#include <complex>
#include "blaze_utils.h"

// Forward Declaration
#include "parameters_fwd.h"

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
    auto psi2 = real(conj(state.psis) % state.psis);
    return psi2 * state.lambda - 1.0;
}

#endif
