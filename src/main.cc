#include <iostream>
#include "common.h"
#include "cosmology.h"
#include "io.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << INFOTAG("Usage: ./sim /path/to/inifile.ini") << std::endl;
        exit(1);
    }

    // Populate Parameter struct
    Parameters param(argv[1]);

    // Initialize the cosmological model, i.e. setup the relation
    // a(tau) and tau(a)
    Cosmology cosmo(param);

    // At this point tau_end is set for the static case as specified in the ini
    // file. This might be wrong, so let's check.
    if (param.cosmo != Parameters::Model::Static)
        param.tau_end = cosmo.tau_of_a(param.a_end);

    // Parameters are now correctly initialized. Print them.
    std::cout << param;

    // Setup HDF5 file
    OutputFile file(param);

#ifndef NDEBUG
    std::cout << INFOTAG("Saving a(tau) to file") << std::endl;
    d_vector tau_grid, a_values;
    double tau = param.tau_start;
    while (tau < param.tau_end) {
        tau_grid.push_back(tau);
        a_values.push_back(cosmo.a_of_tau(tau));
        tau += param.dtau;
    }
    file.write(tau_grid, "tau");
    file.write(a_values, "a_of_tau");
#endif

    // Dummy potential
    SimState state;
    c_vector psis(param.M * param.N, std::complex<double>(1, 1));
    d_vector Vs(param.M * param.N, 0.0);

    state.Vs = Vs;
    state.psis = psis;
    state.n = 1;
    state.tau = 42.2;

    file.write(state, param);
    return 0;
}
