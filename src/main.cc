#include "H5Fpublic.h"
#include "cosmology.h"
#include "hdf5_file.h"
#include "ic.h"
#include "interfaces.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>
#include <chrono>   //runtime measurement
#include <fstream>  //loading the json file
#include <queue>    //priority queue to select to next checkpoint
#include <set>      //building the union of all checkpoints
#include <unordered_map>
#include <vector>

int main(int argc, char **argv) {
    using namespace std::chrono;

#ifdef PSIDMLAB_SMP
    if (!fftw_init_threads()) {
        std::cerr << ERRORTAG("Could not initialize FFTW threading")
                  << std::endl;
        exit(1);
    }
    std::cout << INFOTAG("Threading active") << std::endl;
#endif

    // We need a input file. Otherwise, we give up...
    if (argc != 2) {
        std::cerr << ERRORTAG("Usage: ./sim /path/to/(HDF5|json)") << std::endl;
        exit(1);
    }

    // Parameter datastructure (a json file)
    Parameters param;
    // Load parameters from file
    std::ifstream(argv[1]) >> param;

    const Cosmology cosmo(param);
    param << cosmo;

    // Holds psi, V, n, tau
    SimState state;
    // Populate simulation state
    ICGenerator ic(param);
    ic.generate(state, cosmo);

    const std::string potential = param["Simulation"]["interaction"]["name"];
    auto pot = Interaction::make(potential, param, state);

    const std::string stepper_name = param["Simulation"]["stepper"]["name"];
    auto stepper = TimeEvolution::make(stepper_name, param, state, cosmo);

    // Setup Analysis Functors

    // Mapping from the name of an existing observable to a owning ptr.
    using OwningObservableMap =
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>;
    using ObservableMap = std::unordered_map<std::string, ObservableFunctor *>;

    // A checkpoint consists of its time point, the names of computed
    // observables, and pointers to the latter. A checkpoint does not own
    // the observable !
    using checkpoint_t = std::tuple<double, ObservableMap>;

    // Later on we iterate over checkpoints:
    // For each checkpoint, we...
    // (i)   integrate to the next time point (first tuple element),
    // (ii)  iterate over all observable names (keys in map)
    // (iii) compute the observable via a non-owning pointer (values in map)
    auto cmp = [](const checkpoint_t &a, const checkpoint_t &b) {
        return std::get<0>(a) > std::get<0>(b);
    };
    std::priority_queue<checkpoint_t, std::vector<checkpoint_t>, decltype(cmp)>
        checkpoints(cmp);

    OwningObservableMap observables;
    // Map the name of each observable to its computation time values
    std::unordered_map<std::string, std::vector<double>>
        observable_compute_times;

    // Holds mathematical union of all explicitly stated time points
    std::set<double> time_point_union;
    // Compute time point union and allocate all apearing observables.
    // This populates 'observables' and 'time_point_union'
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
                // Drop all time points that state already visisted (due to ii)
                timepoints.erase(
                    std::remove_if(timepoints.begin(), timepoints.end(),
                                   [&state](const double t) {
                                       return t <= state.tau && state.tau > 0;
                                   }),
                    timepoints.end());

                time_point_union.insert(timepoints.begin(), timepoints.end());
            }
            // Do time points exist after we erased already visited ones?
            if (timepoints.size() != 0) {
                observables[name] =
                    ObservableFunctor::make(ns_name, param, cosmo);
                observable_compute_times[name] = timepoints;
            }
        }
    }
    // Next loop over all unique time points and populate the 'checkpoint' map
    for (double compute_time : time_point_union) {
        ObservableMap obs_at_compute_time;
        for (auto &[name, timepoints] : observable_compute_times) {
            auto obs_ptr = observables[name].get();
            if (timepoints[0] == -1 ||
                std::binary_search(timepoints.begin(), timepoints.end(),
                                   compute_time)) {
                obs_at_compute_time[name] = obs_ptr;
            }
        }
        checkpoints.push(checkpoint_t(compute_time, obs_at_compute_time));
    }

    // Output file setup
    std::string filename = param["General"]["output_file"];

    // File access mode
    HDF5File::Access mode = HDF5File::Access::NewFile;
    if (ic.type == ICType::PreviousSimulation) {
        mode = HDF5File::Access::Write;
    }
    HDF5File file(filename, mode);

    // Define I/O visitor
    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](const auto &output) {
        file.write(path_to_ds, output);
    };

    // Everything is set up...
    // Dump parameters
    std::cout << INFOTAG("Parameter dump:") << std::endl;
    std::cout << std::setw(4) << param << std::endl;

    double runtime_old = 0;

    // For fresh runs (no step performed yet) dump parameters into file
    if (state.tau == 0) {
        file.write_scalar_attribute("/", "parameters", param.dump());
    } else {
        runtime_old = file.read_scalar_attribute<double>("/", "runtime");
    }

    const auto start = high_resolution_clock::now();
    while (!checkpoints.empty()) {
        auto &[checkpoint_time, obs_at_checkpoint] = checkpoints.top();
        stepper->integrate(state, checkpoint_time);
        // Depending on the stepper a retransformation to real space might be
        // necessary. If not, transform acts as identity.
        state.transform(SimState::Representation::Position);
        // Depending on the integrator, state.V might be in an intermediate step
        // and only state.psis is correct. Thus, we recalculate the potential.
        pot->solve(state);

        const double a = cosmo.a_of_tau(checkpoint_time);
        const double z = (a != 0) ? Cosmology::z_of_a(a) : 0;
        std::cout << INFOTAG("Save observables @ tau/z = ")
                  << ((cosmo == CosmoModel::Artificial) ? checkpoint_time : z)
                  << std::flush;
        for (auto [name, obs_ptr] : obs_at_checkpoint) {
            // Construct path to dataset
            auto ds_name = std::to_string(state.n);
            ds_name.insert(ds_name.begin(), 10 - ds_name.length(), '0');
            path_to_ds = "/" + name + "/" + ds_name;

            // Computes the observable variant (i.e. some matrix/vector)
            auto res = obs_ptr->compute(state, observables);
            //            ^                       ^
            //            |       owning ---------|
            //            |---- non-owning
            boost::apply_visitor(write_variant, res);

            // Supplement informations to the observable
            file.write_scalar_attribute(path_to_ds, "tau", state.tau);
            file.write_scalar_attribute(path_to_ds, "a", a);
            file.write_scalar_attribute(path_to_ds, "z", z);
            file.write_scalar_attribute(path_to_ds, "n", state.n);
        }
        std::cout << " ... done" << std::endl;
        checkpoints.pop();
    }

    const auto stop = high_resolution_clock::now();
    const auto runtime = duration_cast<duration<double>>(stop - start);
    file.write_scalar_attribute("/", "runtime", runtime_old + runtime.count());

    return 0;
}
