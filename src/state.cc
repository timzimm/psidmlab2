#include "state.h"
#include "cosmology.h"
#include "parameters.h"

SimState::SimState(const Parameters& p)
    : n(0),
      tau(0),
      dtau(p["Simulation"]["dtau"].get<double>()),
      a(Cosmology::a_of_z(p["Cosmology"]["z_start"].get<double>())) {}

void operator>>(const SimState& state, Parameters& p) {
    p["Simulation"]["N"] = state.psis.rows();
    p["Initial Conditions"]["M"] = state.psis.columns();
}
