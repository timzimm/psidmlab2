#include "cosmology.h"
#include "observables_common.h"

namespace Observable {

// The NLSE of the form
//
//      i∂_t psi = -1/2∂_xx psi + a (G * |psi|^2) psi       (*)
//
// under PERIODIC BOUNDARY CONDITIONS and for symmetric kernels
//
//                  G(x,x') = G(x',x) = G(|x-x'|)
//
// preserves its total linear momentum
//
//                         xmax
//                P(t) = Im ∫dx psi ∂_x psi*                (**)
//                         xmin
//
// The class MomentumDensity computes the NON-INTEGRATED (momentum densities)
// version of (**)
//
//              momentum = Im[ ∂_x psi psi* ]
//
// by a finite difference approximation.

class MomentumDensity : public ObservableFunctor {
    const Domain box;
    double t_prev;
    // Sparse matrix holding the finite difference approximation
    // of the first derivative
    CompressedMatrix<double> grad;
    DynamicVector<double> momentum;

   public:
    MomentumDensity(const Parameters& p, const Cosmology&)
        : box{p}, t_prev(-1), grad(box.N, box.N + 4), momentum(box.N) {
        grad.reserve(4 * box.N);
        // 5 point stencil for first derivative
        for (int i = 0; i < box.N; ++i) {
            grad.append(i, i, 1.0 / (12 * box.dx));
            grad.append(i, i + 1, -2.0 / (3 * box.dx));
            grad.append(i, i + 3, 2.0 / (3 * box.dx));
            grad.append(i, i + 4, -1.0 / (12 * box.dx));
            grad.finalize(i);
        }
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            auto cyclic_extension = [N = box.N](int i) {
                return (N - 2 + i) % N;
            };

            // Cyclic extension of psi to allow gradient computation across
            // boundary
            auto psi = elements(state.psi, cyclic_extension, box.N + 4);
            // Sparse matrix-vector product
            momentum = 0.5 * imag(conj(grad * psi) * state.psi);
        }

        return momentum;
    }
    REGISTER(MomentumDensity)
};

// The class Momentum integrates the momentum density from above
//
//            momentum = Im ∫dx psi ∂_x psi*

class Momentum : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    DynamicVector<double> momentum;

   public:
    Momentum(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_), cosmo(cosmo_), box(p), t_prev(-1), momentum(1, 0) {}

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // Check if there exists a EnergyDensity observable
            const std::string name = "MomentumDensity";
            const std::string full_name = "Observable::" + name;
            auto result = obs.find(name);
            if (result == obs.end()) {
                // Allocate new EnergyDensity observable and
                // store it in global obs map
                obs[name] = ObservableFunctor::make(full_name, p, cosmo);
            }
            ObservableFunctor* pdensity = obs[name].get();
            ReturnType momentum_density = pdensity->compute(state, obs);
            // Trapezodial integration in x-space
            momentum =
                box.dx *
                sum(boost::get<const DynamicVector<double>&>(momentum_density));
        }

        return momentum;
    }
    REGISTER(Momentum)
};

}  // namespace Observable
