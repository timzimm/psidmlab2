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
    // We need a json file. Otherwise, we give up...
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/config_file.json")
                  << std::endl;
        exit(1);
    }

    // Load parameters from file and dump them
    Parameters param;
    std::ifstream(argv[1]) >> param;
    std::cout << INFOTAG("Parsed JSON File. Dump...") << std::endl;
    std::cout << std::setw(4) << param << std::endl;

    // Initialize the cosmological model, i.e. setup the relation
    // a(tau) and tau(a) for the main loop
    Cosmology cosmo(param);

    // Setup initial state
    SimState state(param);
    ICGenerator ic(param);
    ic.generate(state, param);

    // Initialize numerical methods.
    auto integrator_name = param["Simulation"]["integrator"].get<std::string>();
    auto integrator =
        SchroedingerMethod::make(integrator_name, param, state, cosmo);

    // Setup output file
    std::string filename;
    param["General"]["output_file"].get_to(filename);
    HDF5File file(filename);
    file.add_scalar_attribute("/", "parameters", param.dump());

    // Setup Analysis Functors and I/O visitor for their return variant
    std::vector<std::string> keys;
    param["Analysis"]["compute"].get_to(keys);

    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>
        observables;
    for (auto& key : keys)
        observables[key] = ObservableFunctor::make(key, param);

    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](auto& output) {
        file.write(path_to_ds, output);
    };

    // Setup I/O checkpoints
    std::vector<double> save_at;
    param["Analysis"]["save_at"].get_to(save_at);

    // Set final simulation time and checkpoints based on cosmological model
    double tau_end;
    if (cosmo == CosmoModel::Static)
        tau_end = param["Simulation"]["t_end"].get<double>();
    else {
        auto z_end = param["Simulation"]["z_end"].get<double>();
        tau_end = cosmo.tau_of_a(Cosmology::a_of_z(z_end));
        // save_at holds redshifts
        std::transform(std::begin(save_at), std::end(save_at),
                       std::begin(save_at), [&cosmo](double z) {
                           return cosmo.tau_of_a(Cosmology::a_of_z(z));
                       });
    }
    const double dtau = param["Simulation"]["dtau"].get<double>();

    // Everything is set up. Enter the main loop
    auto tau_save = std::begin(save_at);
    bool do_IO = false;

    while (state.tau < tau_end) {
        double tau_dt = state.tau + state.dtau;
        std::cout << "tau_dt " << tau_dt << std::endl;
        std::cout << "tau_save " << *tau_save << std::endl;
        // End of simulation?
        if (tau_dt > tau_end) {
            std::cout << "End" << std::endl;
            state.dtau = tau_end - state.tau;
        }
        // Overstepping I/O checkpoint?
        if (tau_save != std::end(save_at) && tau_dt > *tau_save) {
            std::cout << "I/O" << std::endl;
            state.dtau = *tau_save - state.tau;
            do_IO = true;
        }
        std::cout << "state.tau " << state.tau << std::endl;
        std::cout << "state.dtau " << state.dtau << std::endl;

        integrator->step(state);

        if (do_IO) {
            std::cout << INFOTAG("Save state @ tau = ") << *tau_save
                      << std::endl;
            for (const auto& pair : observables) {
                // Construct path to dataset
                auto& name = pair.first;
                path_to_ds = "/" + name + "/" + std::to_string(state.n);

                // Computes the observable variant;
                auto& routine = pair.second;
                auto res = routine->compute(state);
                boost::apply_visitor(write_variant, res);

                // Supplement informations to the observable
                file.add_scalar_attribute(path_to_ds, "tau", state.tau);
                file.add_scalar_attribute(path_to_ds, "a", state.a);
                file.add_scalar_attribute(path_to_ds, "n", state.n);
            }
            state.dtau = dtau;
            std::next(tau_save);
            do_IO = false;
        }
    }

    return 0;
}
