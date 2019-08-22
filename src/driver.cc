#include "driver.h"
#include "parameters.h"

SimpleDriver::SimpleDriver(const Parameters& p)
    : dt{p["Simulation"]["dtau"].get<double>()} {};

void SimpleDriver::integrate(std::unique_ptr<Stepper>& stepper, SimState& state,
                             const double t_final) {
    while (state.tau + dt < t_final) {
        stepper->step(state, dt);
    }

    stepper->step(state, t_final - state.tau);
}

StableDriver::StableDriver(const Parameters& p){};
void StableDriver::integrate(std::unique_ptr<Stepper>& stepper, SimState& state,
                             const double t_final) {
    for (double dt = stepper->next_dt(state); state.tau + dt < t_final;
         dt = stepper->next_dt(state)) {
        stepper->step(state, dt);
    }

    stepper->step(state, t_final - state.tau);
}
