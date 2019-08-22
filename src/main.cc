#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
#include "parameters.h"
#include "state.h"

#include <chrono>
#include <fstream>
#include <vector>

int main(int argc, char** argv) {
    using namespace std::chrono;

    // We need a json file. Otherwise, we give up...
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/config_file.json")
                  << std::endl;
        exit(1);
    }

    // Load parameters from file
    Parameters param;
    std::ifstream(argv[1]) >> param;

    // Setup a(tau) and tau(a)...
    const Cosmology cosmo(param);
    // ... and transform length scales to code units
    param << cosmo;

    SimState state(param);

    ICGenerator ic(param);
    ic.generate(state, cosmo);
    // Depending on the type of IC chosen, the No. of wavefunctions & No. of
    // dofs (spatial points, fourier modes, basis functions etc.) might have
    // changed. Hence, we inform the parameter file about this potential change.
    state >> param;

    // This defines the PDE to integrate as well as its discretization
    const auto stepper_name = param["Simulation"]["stepper"].get<std::string>();
    auto stepper = Stepper::make(stepper_name, param, state, cosmo);

    // This sets the propagation method, e.g. naive stepping. stability driven,
    // controlled by truncation error etc.
    const auto driver_name = param["Simulation"]["driver"].get<std::string>();
    auto driver = Driver::make(driver_name, param);

    // Setup Analysis Functors
    std::vector<std::string> keys;
    param["Observables"]["compute"].get_to(keys);

    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>
        observables;
    for (auto& key : keys)
        observables[key] = ObservableFunctor::make(key, param, cosmo);

    // Setup I/O checkpoints
    auto save_at = param["Observables"]["save_at"].get<std::vector<double>>();

    // Sort I/O checkpoints and transform them to code time if necessary
    if (cosmo == CosmoModel::Dynamic) {
        auto tau_of_z = [&cosmo](double z) {
            return cosmo.tau_of_a(Cosmology::a_of_z(z));
        };

        // save_at holds redshifts
        std::transform(std::begin(save_at), std::end(save_at),
                       std::begin(save_at), tau_of_z);
    }
    std::sort(save_at.begin(), save_at.end());

    // Output file setup
    std::string filename;
    param["General"]["output_file"].get_to(filename);
    HDF5File file(filename);

    // Define I/O visitor
    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](const auto& output) {
        file.write(path_to_ds, output);
    };

    // Everything is set up...
    // Dump parameters (modulo the i/o checkpoints since this can be a long
    // list)
    param["Observables"]["save_at"] = "";
    std::cout << INFOTAG("Parameter dump:") << std::endl;
    std::cout << std::setw(4) << param << std::endl;
    file.add_scalar_attribute("/", "parameters", param.dump());

    const auto start = high_resolution_clock::now();
    for (auto checkpoint : save_at) {
        driver->integrate(stepper, state, checkpoint);
        // Depending on the stepper a retransformation to real space might be
        // necessary. If not, transform acts as identity.
        state.transform(SimState::Representation::Position);

        const double a = cosmo.a_of_tau(checkpoint);
        const double z = Cosmology::z_of_a(a);
        std::cout << INFOTAG("Save state @ tau/z = ")
                  << ((cosmo == CosmoModel::Static) ? checkpoint : z)
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
            file.add_scalar_attribute(path_to_ds, "a", a);
            file.add_scalar_attribute(path_to_ds, "z", z);
            file.add_scalar_attribute(path_to_ds, "n", state.n);
        }
        std::cout << " ... done" << std::endl;
    }

    const auto stop = high_resolution_clock::now();
    const auto runtime = duration_cast<duration<double>>(stop - start);
    file.add_scalar_attribute("/", "runtime", runtime.count());

    return 0;
}
