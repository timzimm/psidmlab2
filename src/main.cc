#include <blaze/math/Row.h>
#include <iostream>
#include <memory>
#include "common.h"
#include "cosmology.h"
#include "fftw_alloc.h"
#include "ic.h"
#include "io.h"
#include "poisson_solver.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/inifile.ini") << std::endl;
        exit(1);
    }

    // Populate Parameter struct
    Parameters param(argv[1]);

    // Initialize the cosmological model, i.e. setup the relation
    // a(tau) and tau(a)
    Cosmology cosmo(param);

    // At this point tau_end is set for the static case as specified in the
    // ini file. This might be wrong, so let's check.
    if (param.cosmo != Parameters::Model::Static)
        param.tau_end = cosmo.tau_of_a(param.a_end);

    // Parameters are now correctly initialized. Print them.
    std::cout << param;

    SimState state(param);

    ICGenerator ic(param);
    ic.generate(state);

    // Setup HDF5 file
    OutputFile file(param);

#ifndef NDEBUG
    /* file.write(blaze::trans(blaze::row(state.Vs, 0)), "rho"); */

    // Poisson Solver Test
    param.N = 256;
    param.L = 2 * M_PI;
    param.dx = param.L / param.N;
    blaze::DynamicVector<double, blaze::rowVector> source(param.N);

    std::unique_ptr<Solvers::Poisson> poisson =
        std::make_unique<Solvers::FFT>(param);
    for (int j = 0; j < 10; ++j) {
        double delta = param.L / 20 * j;
        for (int i = 0; i < source.size(); ++i)
            source[i] = delta + param.dx * i;
        source = blaze::cos(source);
        blaze::subvector(source, 10, 20) = 1;
        blaze::row(state.Vs, 2 * j) = source;
        blaze::row(state.Vs, 2 * j + 1) = poisson->solve(source);
        std::string sou = std::string("source") + std::to_string(j);
        std::string rho = std::string("rho") + std::to_string(j);
        file.write(blaze::trans(blaze::row(state.Vs, 2 * j)), sou);
        file.write(blaze::trans(blaze::row(state.Vs, 2 * j + 1)), rho);
    }

    /* file.write(state.lambda, "lambda"); */
    /* std::cout << INFOTAG("Saving a(tau) to file") << std::endl; */
    /* std::vector<double> tau_grid, a_values; */
    /* double tau = param.tau_start; */
    /* while (tau < param.tau_end) { */
    /*     tau_grid.push_back(tau); */
    /*     a_values.push_back(cosmo.a_of_tau(tau)); */
    /*     tau += param.dtau; */
    /* } */
    /* file.write(tau_grid, "tau"); */
    /* file.write(a_values, "a_of_tau"); */
#endif

    return 0;
}
