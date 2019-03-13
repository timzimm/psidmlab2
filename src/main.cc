#include <iostream>
#include "common.h"
#include "io.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << INFOTAG("Usage: ./sim /path/to/inifile.ini") << std::endl;
        exit(1);
    }

    Parameters parameters(argv[1]);
    std::cout << parameters;

    // Dummy potential
    SimState state;
    c_vector psis(parameters.M * parameters.N, std::complex<double>(1, 1));
    d_vector Vs(parameters.M * parameters.N, 0.0);

    state.Vs = Vs;
    state.psis = psis;
    state.n = 1;
    state.tau = 42.2;

    OutputFile file_hdf5(parameters);
    file_hdf5.write(state, parameters);
    return 0;
}
