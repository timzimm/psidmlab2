#include <iostream>
#include <memory>
#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

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
    /* if (param.cosmo != CosmoModel::Static) */
    /* param.tau_end = cosmo.tau_of_a(param.a_end); */

    // Parameters are now correctly initialized. Print them.
    std::cout << param;

    SimState state(param);

    ICGenerator ic(param);
    ic.generate(state);

    // Setup HDF5 file
    OutputFile file(param);

    // Initiliaze numerical method. This selects both the Schroedinger and
    // Poisson algorithm to be used for the integration
    /* auto integrator = SchroedingerMethod::make(param.integrator, param); */

    /* while (state.tau < tau_end) { */
    /*     (*integrator)(state); */
    /* } */

    /* auto psi = blaze::row(state.psis, 1); */
    /* auto psi_source = blaze::real(psi); */
    /* file.write(blaze::trans(psi_source), "psi2"); */

    file.write(state.V, "V");
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

    return 0;
}
