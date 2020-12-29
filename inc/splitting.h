#ifndef __SPLITTING__
#define __SPLITTING__

#include "driver.h"
#include "parameters.h"

template <typename FlowMapA, typename FlowMapB>
class PRK : public DefaultDriver<PRK<FlowMapA, FlowMapB>> {
    // Subproblem time evolution operators (A, B possible nonlinear)
    FlowMapA phiA; // integrator for del_t psi = A(psi)
    FlowMapB phiB; // integrator for del_t psi = B(psi)
    const int order;
    blaze::DynamicVector<double> a;
    blaze::DynamicVector<double> b;
    const bool stable;
    double dt;

  public:
    PRK(const Parameters &p, const SimState &state, const Cosmology &cosmo_)
        : DefaultDriver<PRK<FlowMapA, FlowMapB>>(p), phiA(p, state, cosmo_),
          phiB(p, state, cosmo_), order{p["Simulation"]["stepper"]["order"]},
          stable{p["Simulation"]["driver"]["stable"]},
          dt{p["Simulation"]["driver"]["dtau"].get<double>()} {
        switch (order) {
        case 2:
            // Strang Splitting
            b = {0.5, 0.5};
            a = {1};
            break;
        case 4:
            // Blanes, Moan symmetric PRK splitting
            a = {
                0.209515106613361891,  -0.143851773179818077,
                0.434336666566456186,  0.434336666566456186,
                -0.143851773179818077, 0.209515106613361891,
            };
            b = {
                0.0792036964311954608,  0.353172906049773948,
                -0.0420650803577191948, 0.219376955753499572,
                -0.0420650803577191948, 0.353172906049773948,
                0.0792036964311954608,
            };
            break;
        default:
            std::cout << ERRORTAG("Order not supported") << std::endl;
            exit(1);
        }
    }
    // Next time step for total problem determined by stability of the
    // subproblems A, B and the artificial upper bound dt
    double next_dt(const SimState &state) const {
        return std::min({phiA.next_dt(state), phiB.next_dt(state), dt});
    }

    // Evolution operator for i del_t psi = (A + B) psi
    void step(SimState &state, const double dt) const {
        const double t = state.tau;
        const int n = state.n;
        for (int stage = 0; stage < a.size(); ++stage) {
            phiB.step(state, b[stage] * dt);
            state.tau = t;
            phiA.step(state, a[stage] * dt);
            state.tau = t;
        }

        phiB.step(state, b[b.size() - 1] * dt);
        state.tau = t + dt;
        state.n = n + 1;
    }

    /*
    void integrate(SimState& state, const double t_final) {
        if (stable) {
            for (double dt = next_dt(state); state.tau + dt < t_final;
                 dt = next_dt(state)) {
                step(state, dt);
            }
        } else {
            double t = state.tau;
            int n = state.n;
            if (state.tau + 2 * dt < t_final) {
                // pre-step - omit last applitcation B flow
                for (int stage = 0; stage < a.size(); ++stage) {
                    phiB.step(state, b[stage] * dt);
                    state.tau = t;
                    state.n = n;
                    phiA.step(state, a[stage] * dt);
                    state.tau = t;
                    state.n = n;
                }
                while (state.tau + dt < t_final) {
                    // Join adjacent B-flows across steps (FSAL)
                    phiB.step(state, 2 * b[0] * dt);
                    state.tau = t;
                    state.n = n;
                    phiA.step(state, a[0] * dt);
                    state.tau = t;
                    state.n = n;
                    for (int stage = 1; stage < a.size(); ++stage) {
                        phiB.step(state, b[stage] * dt);
                        state.tau = t;
                        state.n = n;
                        phiA.step(state, a[stage] * dt);
                        state.tau = t;
                        state.n = n;
                    }
                    t = state.tau = t + dt;
                    n = state.n += 1;
                }
                // post step - synchonize with time grid by carrying
                // out the residual step
                phiB.step(state, b[b.size() - 1] * dt);
                state.tau = t;
                state.n = n;
            } else if (state.tau + dt < t_final) {
                step(state, dt);
            }

            // Residual step to requested final time
            step(state, t_final - state.tau);
        }
    }
    */
    REGISTER(PRK)
};

#endif
