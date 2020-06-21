#include "evolution/interaction_potential_cap_magnus.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

#include <fstream>

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotentialCAPMagnus::InteractionPotentialCAPMagnus(
    const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo(cosmo_),
      box(p),
      dt_last(-1),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)},
      phi(box.N),
      attenuator(box.N) {
    std::ifstream cap_file(
        p["Simulation"]["stepper"]["cap_file"].get<std::string>());
    int data_N = std::count(std::istreambuf_iterator<char>(cap_file),
                            std::istreambuf_iterator<char>(), '\n');
    cap_file.seekg(0);
    if (data_N != box.N) {
        std::cerr << ERRORTAG("#lines in cap_file(" << data_N << ") != N")
                  << std::endl;
        exit(1);
    }
    fill_from_file(cap_file, phi);
}

// Non linear phase method
double InteractionPotentialCAPMagnus::next_dt(const SimState& state) const {
    const double phi_max = 0.1;
    return phi_max / (cosmo.a_of_tau(state.tau) * max(abs(state.V)));
}

void InteractionPotentialCAPMagnus::step(SimState& state,
                                         const double dt) const {
    state.transform(SimState::Representation::Position);

    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    if (dt_last - dt != 0) {
        attenuator = exp(-phi * dt / 2.0);
    }

    // In place computation of potential
    state.V = attenuator * rho_from(state);
    pot->solve(state.V, state.V);

    state.psi *= exp(-1.0i * a_half * state.V * dt) * attenuator;

    state.tau += dt;
    state.n += 1;
    dt_last = dt;
}

}  // namespace Schroedinger
