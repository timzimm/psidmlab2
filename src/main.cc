#include <iostream>
#include "io.h"
#include "psidm.h"

int main() {
    Parameters parameters;
    parameters.filename = "sim.h5";
    parameters.M = 5;
    parameters.N = 10;

    // Dummy potential
    SimState state;
    c_vector psis(parameters.M * parameters.N, std::complex<double>(1, 1));
    d_vector Vs(parameters.M * parameters.N, 0.0);
    /* for (int i = 0; i < parameters.M; ++i) { */
    /*     for (int j = 0; j < parameters.N; ++j) { */
    /*         Vs[i][j] += i; */
    /*         /1* psis[i][j] += std::complex<double>(0, i); *1/ */
    /*         /1* std::cout << Vs[i][j] << std::endl; *1/ */
    /*     } */
    /* } */

    state.Vs = Vs;
    state.psis = psis;
    state.n = 1;
    state.t = 42.2;

    OutputFile file_hdf5(parameters);
    file_hdf5.write(state, parameters);
    return 0;
}
