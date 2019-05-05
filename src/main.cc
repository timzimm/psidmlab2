#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>
#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/config_file.json")
                  << std::endl;
        exit(1);
    }

    Parameters param;
    std::ifstream(argv[1]) >> param;
    std::cout << INFOTAG("Parsed JSON File. Dump...") << std::endl;
    std::cout << std::setw(4) << param << std::endl;

    // Initialize the cosmological model, i.e. setup the relation
    // a(tau) and tau(a) for the main loop
    Cosmology cosmo(param);

    // At this point tau_end is set for the static case as specified in the
    // ini file. This might be wrong, so let's check.
    /* if (param.cosmo != CosmoModel::Static) */
    /* param.tau_end = cosmo.tau_of_a(param.a_end); */

    // Setup initial state
    SimState state(param);
    ICGenerator ic(param);
    ic.generate(state);

    // Setup Analysis Functors
    std::string filename;
    param["General"]["output_file"].get_to(filename);

    std::vector<std::string> keys;
    param["Analysis"]["compute"].get_to(keys);

    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>
        observables;
    for (auto& key : keys)
        observables[key] = ObservableFunctor::make(key, param);

    HDF5File file(filename);

    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](auto& output) {
        file.write(path_to_ds, output);
    };

    std::string parameter_dump = param.dump();

    file.add_scalar_attribute("/", "parameters", parameter_dump);

    for (const auto& pair : observables) {
        // Construct path to dataset
        auto& name = pair.first;
        path_to_ds = "/" + name + "/0";

        // Computes the observable variant;
        auto& routine = pair.second;
        auto res = routine->compute(state);
        boost::apply_visitor(write_variant, res);

        // Supplement informations to the observable
        std::map<std::string, double> m{{"t", 14.5}};
        file.add_scalar_attribute(path_to_ds, m);
    }

    // Initiliaze numerical method. This selects both the Schroedinger and
    // Poisson algorithm to be used for the integration
    /* auto integrator = SchroedingerMethod::make(param.integrator, param); */

    /* while (state.tau < tau_end) { */
    /*     (*integrator)(state); */
    /* } */

    /* auto psi = blaze::row(state.psis, 1); */
    /* auto psi_source = blaze::real(psi); */
    /* file.write(blaze::trans(psi_source), "psi2"); */

    /* file.write(state.V, "V"); */
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
