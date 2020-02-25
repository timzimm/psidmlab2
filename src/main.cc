#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
#include "parameters.h"
#include "state.h"

#include <chrono>   //runtime measurement
#include <fstream>  //loading the json file
#include <queue>    //priority queue to select to next checkpoint
#include <set>      //building the union of all checkpoints
#include <vector>

int main(int argc, char **argv) {
    using namespace std::chrono;

#ifdef PSIMDLAB_SMP
    if (!fftw_init_threads()) {
        std::cerr << ERRORTAG("Could not initialize FFTW threading")
                  << std::endl;
        exit(1);
    }
#endif

    // We need a json file. Otherwise, we give up...
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/config_file.json")
                  << std::endl;
        exit(1);
    }

    // Load parameters from file
    Parameters param;
    std::ifstream(argv[1]) >> param;

    // Holds mathematical union of all explicitly stated time points
    std::set<double> time_point_union;
    for (auto &[name, parameters] : param["Observables"].items()) {
        if (auto compute_at = parameters["compute_at"];
            compute_at != Parameters::array() &&
            compute_at != Parameters::array({-1})) {
            // Extract time points at which the observable should be saved
            auto timepoints = compute_at.get<std::vector<double>>();
            time_point_union.insert(timepoints.begin(), timepoints.end());
        }
    }
    param["Cosmology"]["save_at"] = time_point_union;
    time_point_union.clear();

    // Setup a(tau) and tau(a)...
    const Cosmology cosmo(param["Cosmology"]);
    // ... and transform length scales to code units
    param << cosmo;

    SimState state(param);

    ICGenerator ic(param);
    ic.generate(state, cosmo);
    state >> param;

    const auto potential = param["Simulation"]["potential"].get<std::string>();
    auto pot = Interaction::make(potential, param, state);

    // This defines the PDE to integrate as well as its discretization
    const auto stepper_name =
        param["Simulation"]["stepper"]["name"].get<std::string>();
    auto stepper = TimeEvolution::make(stepper_name, param, state, cosmo);

    // Setup Analysis Functors
    using observable_ptr = std::unique_ptr<ObservableFunctor>;
    // A checkpoint consists of its time point , the names of computed
    // observables, and pointers to the latter
    using checkpoint_t = std::tuple<double, std::vector<std::string>,
                                    std::vector<observable_ptr *>>;

    std::unordered_map<std::string,
                       std::pair<std::vector<double>, observable_ptr>>
        observables;
    std::priority_queue<checkpoint_t, std::vector<checkpoint_t>,
                        std::greater<checkpoint_t>>
        checkpoints;

    for (auto &[name, parameters] : param["Observables"].items()) {
        auto compute_at = parameters["compute_at"];
        if (compute_at != Parameters::array()) {
            // Full qualified name of observable class
            const std::string ns_name = "Observable::" + name;
            // Extract time points at which the observable should be saved
            auto timepoints = compute_at.get<std::vector<double>>();
            // only explicitly stated checkpoints contribute to the union
            if (compute_at != Parameters::array({-1})) {
                // transform to code time if necessary
                if (cosmo == CosmoModel::Dynamic) {
                    auto tau_of_z = [&cosmo](double z) {
                        return cosmo.tau_of_a(Cosmology::a_of_z(z));
                    };
                    // save_at holds redshifts
                    std::transform(std::begin(timepoints), std::end(timepoints),
                                   std::begin(timepoints), tau_of_z);
                }
                time_point_union.insert(timepoints.begin(), timepoints.end());
            }
            observables[name] =
                std::make_pair(std::move(timepoints),
                               ObservableFunctor::make(ns_name, param, cosmo));
        }
        compute_at = "";
    }
    for (double checkpoint : time_point_union) {
        std::vector<std::string> obs_names;
        std::vector<observable_ptr *> obs_ptrs;
        for (auto &[name, value] : observables) {
            auto &[timepoints, obs_ptr] = value;
            if (timepoints[0] == -1 ||
                std::binary_search(timepoints.begin(), timepoints.end(),
                                   checkpoint)) {
                obs_ptrs.push_back(&obs_ptr);
                obs_names.push_back(name);
            }
        }
        checkpoints.push(checkpoint_t(checkpoint, obs_names, obs_ptrs));
    }

    // Output file setup
    std::string filename;
    param["General"]["output_file"].get_to(filename);
    HDF5File file(filename);

    // Define I/O visitor
    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](const auto &output) {
        file.write(path_to_ds, output);
    };

    // Everything is set up...
    // Dump parameters (modulo the i/o checkpoints since this can be a long
    // list)
    std::cout << INFOTAG("Parameter dump:") << std::endl;
    std::cout << std::setw(4) << param << std::endl;
    file.add_scalar_attribute("/", "parameters", param.dump());

    const auto start = high_resolution_clock::now();
    while (!checkpoints.empty()) {
        auto &[checkpoint, names, obs_ptrs] = checkpoints.top();
        stepper->integrate(state, checkpoint);
        // Depending on the stepper a retransformation to real space might be
        // necessary. If not, transform acts as identity.
        state.transform(SimState::Representation::Position);
        // Depending on the integrator, state.V might be in an intermediate step
        // and only state.psis is correct. Thus, we recalculate the potential.
        pot->solve(state);

        const double a = cosmo.a_of_tau(checkpoint);
        const double z = Cosmology::z_of_a(a);
        std::cout << INFOTAG("Save observables @ tau/z = ")
                  << ((cosmo == CosmoModel::Artificial) ? checkpoint : z)
                  << std::flush;
        for (int i = 0; i < obs_ptrs.size(); ++i) {
            auto &obs_ptr = obs_ptrs[i];
            auto &name = names[i];
            // Construct path to dataset
            auto ds_name = std::to_string(state.n);
            ds_name.insert(ds_name.begin(), 10 - ds_name.length(), '0');
            path_to_ds = "/" + name + "/" + ds_name;

            // Computes the observable variant (i.e. some matrix/vector)
            auto res = (*obs_ptr)->compute(state);
            boost::apply_visitor(write_variant, res);

            // Supplement informations to the observable
            file.add_scalar_attribute(path_to_ds, "tau", state.tau);
            file.add_scalar_attribute(path_to_ds, "a", a);
            file.add_scalar_attribute(path_to_ds, "z", z);
            file.add_scalar_attribute(path_to_ds, "n", state.n);
        }
        std::cout << " ... done" << std::endl;
        checkpoints.pop();
    }

    const auto stop = high_resolution_clock::now();
    const auto runtime = duration_cast<duration<double>>(stop - start);
    file.add_scalar_attribute("/", "runtime", runtime.count());

    return 0;
}
