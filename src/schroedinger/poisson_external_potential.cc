#include "schroedinger/poisson_external_potential.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

#include <fstream>

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

PoissonExternalPotential::PoissonExternalPotential(const Parameters& p,
                                                   const SimState& state,
                                                   const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      pot_external(0) {
    const int N = p["Simulation"]["N"].get<int>();
    std::ifstream pot_file{
        p["Simulation"]["stepper"]["external_potential"].get<std::string>()};
    if (!pot_file) {
        std::cout << ERRORTAG("external potential file not found") << std::endl;
        exit(1);
    }

    // Determine # rows in file
    double dataN = std::count(std::istreambuf_iterator<char>(pot_file),
                              std::istreambuf_iterator<char>(), '\n');
    pot_file.seekg(0);
    if (dataN != N) {
        std::cout << ERRORTAG("N from external potential differs") << std::endl;
        exit(1);
    }
    pot_external.resize(N);
    fill_from_file(pot_file, pot_external);
}

// Non linear phase method with phi_max = pi/2
double PoissonExternalPotential::next_dt(const SimState& state) const {
    return M_PI / (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void PoissonExternalPotential::step(SimState& state, const double dt) {
    state.transform(SimState::Representation::Position);
    // |psi|^2 is invariant. Thus we can use the input state to determine the
    // gravitational potential for the entire time step.
    pot->solve(state);
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    auto drift = expand(
        exp(-1.0i * (a_half * state.V * dt + pot_external * dt)), state.M);
    state.psis = drift % state.psis;

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
