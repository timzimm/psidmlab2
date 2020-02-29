#include "H5Fpublic.h"
#include "cosmology.h"
#include "ic.h"
#include "interfaces.h"
#include "io.h"
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

    // Two operationn modi need to be distinguish:
    // (i) Fresh simulation
    // (ii) Resumed simulation
    bool resume = false;

    // Holds psi, V, n, tau
    SimState state;

    // Parameter datastructure (a json file)
    Parameters param;

    // Get simulation parameters via input file
    // The input file is either ...
    // (i) json file for a fresh simulation.
    // (ii) a HDF5 file holding an existing simulation state to be resumed
    auto result = H5Fis_hdf5(argv[1]);
    if (result < 0) {
        std::cerr
            << ERRORTAG(
                   "Cannot determine input file type. Does the file exist?")
            << std::endl;
        // Case (i) :  New simulation
    } else if (result == 0) {
        std::cout << "Initialize simulation from json file" << std::endl;
        // Load parameters from file
        std::ifstream(argv[1]) >> param;
        // Case (ii) :  Resume old simulation
    } else {
        resume = true;
        std::cout << "Resume simulation with last state stored in input file"
                  << std::endl;

        HDF5File file(argv[1]);
        // Load parameters from dumped json in HDF5 file
        std::string dump = file.read_string_attribute("/", "parameters");
        param = Parameters::parse(dump);

        // Determine last available wavefunction
        std::vector<std::string> psis = file.ls("/WaveFunction");
        auto psi0 = std::max_element(psis.begin(), psis.end());
        if (psi0 != psis.end()) {
            auto psi0_matrix = file.read_matrix(*psi0);
            // Populate simualation state.
            // We don't need the potential. It can be computed from psi.
            state.psi = blaze::map(blaze::column(psi0_matrix, 0),
                                   blaze::column(psi0_matrix, 1),
                                   [](const double r, const double i) {
                                       return std::complex<double>(r, i);
                                   });
            state.tau = std::stod(file.read_string_attribute(*psi0, "tau"));
            state.n = std::stoull(file.read_string_attribute(*psi0, "n"));
        } else {  // If there is no wave function we're back to (i)
            resume = false;
        }
    }

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

    const Cosmology cosmo(param["Cosmology"]);
    param << cosmo;

    // In case (i) initial conditions are generated according to config file.
    if (!resume) {
        ICGenerator ic(param);
        // Populates simulation state
        ic.generate(state, cosmo);
        state >> param;
    }

    const auto potential = param["Simulation"]["potential"].get<std::string>();
    auto pot = Interaction::make(potential, param, state);

    // This defines the PDE to integrate as well as its discretization
    const auto stepper_name =
        param["Simulation"]["stepper"]["name"].get<std::string>();
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
                    std::remove_if(
                        timepoints.begin(), timepoints.end(),
                        [&state](const double t) { return t < state.tau; }),
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
    std::string filename;
    param["General"]["output_file"].get_to(filename);
    HDF5File file(filename);

    // Define I/O visitor
    std::string path_to_ds;
    auto write_variant = [&file, &path_to_ds](const auto &output) {
        file.write(path_to_ds, output);
    };

    // Everything is set up...
    // Dump parameters
    std::cout << INFOTAG("Parameter dump:") << std::endl;
    std::cout << std::setw(4) << param << std::endl;
    if (!resume) {
        file.add_scalar_attribute("/", "parameters", param.dump());
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
        const double z = Cosmology::z_of_a(a);
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
