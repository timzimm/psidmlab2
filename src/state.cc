#include "state.h"
#include "cosmology.h"
#include "parameters.h"

SimState::SimState(const Parameters& p)
    : n(0),
      tau(0),
      dtau(p.get<double>("dtau")),
      a(Cosmology::a_of_z(p.get<double>("z_start"))) {}
