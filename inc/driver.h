#ifndef __DRIVER__
#define __DRIVER__
#include "interfaces.h"
#include "parameters.h"

template <typename Derived>
class DefaultDriver : public TimeEvolution {
    const bool stable;
    double dt;

   public:
    DefaultDriver(const Parameters& p)
        : stable{p["Simulation"]["driver"]["stable"]},
          dt{stable ? 0 : p["Simulation"]["driver"]["dtau"].get<double>()} {}

    const Derived& self() const { return static_cast<const Derived&>(*this); }
    Derived& self() { return static_cast<Derived&>(*this); }

    // Forward to true implementation
    void step(SimState& state, const double dt) const {
        self().step(state, dt);
    }
    double next_dt(const SimState& state) { return self().next_dt(state); }
    // Default implementation of integration functionality. This is used if
    // Derived does not explicitly overload integrate(...)
    void integrate(SimState& state, const double t_final) {
        if (stable) {
            for (double dt = next_dt(state); state.tau + dt < t_final;
                 dt = next_dt(state)) {
                step(state, dt);
            }
        } else {
            while (state.tau + dt < t_final) {
                step(state, dt);
            }
        }
        // Residual step to exactly arrive at the final time
        step(state, t_final - state.tau);
    }
};

#endif
