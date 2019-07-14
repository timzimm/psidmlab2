#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
#include "logging.h"
#include "parameters.h"
#include "state.h"

#include <chrono>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <vector>

int main(int argc, char** argv) {
    using namespace std::chrono;

    // We need a json file. Otherwise, we give up...
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/config_file.json")
                  << std::endl;
        exit(1);
    }

    // Load parameters from file and dump them
    Parameters param;
    std::ifstream(argv[1]) >> param;

    // Initialize the cosmological model, i.e. setup the relation
    // a(tau) and tau(a) for the main loop
    const Cosmology cosmo(param);

    SimState state(param);
    ICGenerator ic(param);
    ic.generate(state, cosmo);

    const auto integrator_name =
        param["Simulation"]["integrator"].get<std::string>();
    auto integrator =
        SchroedingerMethod::make(integrator_name, param, state, cosmo);

    // Bring parameter datastructure back to a consistent state (N, M)
    state >> param;
    std::cout << INFOTAG("Parameter dump:") << std::endl;
    std::cout << std::setw(4) << param << std::endl;

    std::string filename;
    param["General"]["output_file"].get_to(filename);
    HDF5File file(filename);
    file.add_scalar_attribute("/", "parameters", param.dump());

    // Setup Analysis Functors and I/O visitor for their return variant
    std::vector<std::string> keys;
    param["Observables"]["compute"].get_to(keys);

    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>
        observables;
    for (auto& key : keys)
        observables[key] = ObservableFunctor::make(key, param);

    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](const auto& output) {
        file.write(path_to_ds, output);
    };

    // Setup I/O checkpoints
    auto save_at = param["Observables"]["save_at"].get<std::vector<double>>();

    // Set final simulation time and checkpoints based on cosmological model
    const double tau_end = [&] {
        if (cosmo == CosmoModel::Dynamic) {
            auto tau_of_a = [&cosmo](double z) {
                return cosmo.tau_of_a(Cosmology::a_of_z(z));
            };

            // save_at holds redshifts
            std::transform(std::begin(save_at), std::end(save_at),
                           std::begin(save_at), tau_of_a);
        }
        return *std::max_element(std::begin(save_at), std::end(save_at));
    }();

    const double dtau = param["Simulation"]["dtau"].get<double>();

    auto tau_save = std::begin(save_at);
    bool do_IO = false;

    const auto start = high_resolution_clock::now();

    while (state.tau <= tau_end) {
        double tau_dt = state.tau + state.dtau;

        // Overstepping I/O checkpoint?
        if (tau_save != std::end(save_at) && tau_dt >= *tau_save) {
            state.dtau = *tau_save - state.tau;
            do_IO = true;
        }

        integrator->step(state);

        if (do_IO) {
            std::cout << INFOTAG("Save state @ tau/z=")
                      << ((cosmo == CosmoModel::Static)
                              ? *tau_save
                              : Cosmology::z_of_a(cosmo.a_of_tau(*tau_save)))
                      << std::flush;
            for (const auto& pair : observables) {
                // Construct path to dataset
                auto& name = pair.first;
                auto ds_name = std::to_string(state.n);
                ds_name.insert(ds_name.begin(), 10 - ds_name.length(), '0');
                path_to_ds = "/" + name + "/" + ds_name;

                // Computes the observable variant (i.e. some matrix/vector)
                auto& routine = pair.second;
                auto res = routine->compute(state);
                boost::apply_visitor(write_variant, res);

                // Supplement informations to the observable
                file.add_scalar_attribute(path_to_ds, "tau", state.tau);
                file.add_scalar_attribute(path_to_ds, "a", state.a);
                file.add_scalar_attribute(path_to_ds, "z",
                                          Cosmology::z_of_a(state.a));
                file.add_scalar_attribute(path_to_ds, "n", state.n);
            }
            std::cout << " ... done" << std::endl;
            state.dtau = dtau;
            tau_save++;
            do_IO = false;
        }
    }

    const auto stop = high_resolution_clock::now();
    const auto runtime = duration_cast<duration<double>>(stop - start);
    file.add_scalar_attribute("/", "runtime", runtime.count());

    return 0;
}
