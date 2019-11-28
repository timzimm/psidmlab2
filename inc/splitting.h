#ifndef __SPLITTING__
#define __SPLITTING__

#include "driver.h"
#include "parameters.h"

template <typename FlowMapA, typename FlowMapB>
class SRKN : public DefaultDriver<SRKN<FlowMapA, FlowMapB>> {
    // Subproblem Integrators (A, B possible nonlinear)
    FlowMapA phiA;  // integrator for del_t psi = A psi
    FlowMapB phiB;  // integrator for del_t psi = B psi
    int order;
    blaze::DynamicVector<double> a;
    blaze::DynamicVector<double> b;
    bool stable;
    double dt;

   public:
    SRKN(const Parameters& p, const SimState& state, const Cosmology& cosmo_)
        : DefaultDriver<SRKN<FlowMapA, FlowMapB>>(p),
          phiA(p, state, cosmo_),
          phiB(p, state, cosmo_),
          order{p["Simulation"]["stepper"]["order"].get<int>()},
          stable{p["Simulation"]["driver"]["stable"].get<bool>()},
          dt{p["Simulation"]["driver"]["dtau"].get<double>()} {
        switch (order) {
            case 2:
                // Strang Splitting
                b = {0.5, 0.5};
                a = {1};
                break;
            case 4:
                // Blanes, Moan SRKN_6 BAB splitting
                b = {.0829844064174052,
                     .39630980149836,
                     -.0390563049223486,
                     0,
                     0,
                     0,
                     0};
                b[3] = 1 - 2 * (b[0] + b[1] + b[2]);
                subvector(b, 4, 3) = reverse(subvector(b, 0, 3));

                a = {0.2452989571842710, 0.6048726657110800, 0, 0, 0, 0};
                a[2] = 0.5 - (a[0] + a[1]);
                subvector(a, 3, 3) = reverse(subvector(a, 0, 3));
                break;
            case 6:
                // Blanes, Moan SRKN_11 BAB splitting
                b = {.0414649985182624,
                     .198128671918067,
                     -.0400061921041533,
                     .0752539843015807,
                     -.0115113874206879,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0};
                b[5] = 0.5 - (b[0] + b[1] + b[2] + b[3] + b[4]);
                subvector(b, 5, 6) = reverse(subvector(b, 0, 6));

                a = {.123229775946271,
                     .290553797799558,
                     -.127049212625417,
                     -.246331761062075,
                     .357208872795928,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0};
                a[5] = 1 - 2 * (a[0] + a[1] + a[2] + a[3] + a[4]);
                subvector(a, 6, 5) = reverse(subvector(b, 0, 5));
                break;
            default:
                std::cout << ERRORTAG("Order not supported") << std::endl;
                exit(1);
        }
    }
    // Next time step for total problem determined by stability of the
    // subproblems A and B, i.e. taking the minium of both
    double next_dt(const SimState& state) const {
        return std::min({phiA.next_dt(state), phiB.next_dt(state)});
    }

    // Evolution operator for i del_t psi = (A + B) psi
    void step(SimState& state, const double dt) {
        const double t = state.tau;
        const int n = state.n;
        for (int stage = 0; stage < a.size(); ++stage) {
            phiB.step(state, b[stage] * dt);
            state.tau = t;
            state.n = n;
            phiA.step(state, a[stage] * dt);
            state.tau = t;
            state.n = n;
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
    REGISTER(SRKN)
};

#endif
