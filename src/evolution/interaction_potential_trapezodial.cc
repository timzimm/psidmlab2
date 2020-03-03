#include "evolution/interaction_potential_trapezodial.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

#include <fstream>

namespace Schroedinger {

using namespace blaze;
using namespace std::complex_literals;

InteractionPotentialTrapezodial::InteractionPotentialTrapezodial(
    const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo{cosmo_},
      N{p["Simulation"]["N"]},
      pot{Interaction::make(p["Simulation"]["potential"], p, state)},
      pot_external(0) {
    try {
        std::ifstream pot_file{p["Simulation"]["stepper"]["external_potential"]
                                   .get<std::string>()};

        if (!pot_file) return;
        std::cout << INFOTAG("External potential loaded") << std::endl;

        // Determine # rows in file
        double dataN = std::count(std::istreambuf_iterator<char>(pot_file),
                                  std::istreambuf_iterator<char>(), '\n');
        pot_file.seekg(0);
        if (dataN != N) {
            std::cout << ERRORTAG("N from external potential differs")
                      << std::endl;
            exit(1);
        }
        pot_external.resize(N);
        fill_from_file(pot_file, pot_external);
    } catch (nlohmann::detail::type_error) {
    }
}

// Non linear phase method with phi_max = pi/2
double InteractionPotentialTrapezodial::next_dt(const SimState& state) const {
    return M_PI /
           (2 * cosmo.a_of_tau(state.tau) * max(abs(state.V + pot_external)));
}

void InteractionPotentialTrapezodial::step(SimState& state,
                                           const double dt) const {
    state.transform(SimState::Representation::Position);
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);
    pot->solve(state);
    if (pot_external.size() == N) state.V += pot_external;
    state.psi *= exp(-1.0i * a_half * state.V * dt);

    state.tau += dt;
    state.n += 1;
}

}  // namespace Schroedinger
