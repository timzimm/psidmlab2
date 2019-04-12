#ifndef __STATE__
#define __STATE__
#include <complex>
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"

// Forward Declaration
class Parameters;

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

#endif
