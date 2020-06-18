#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"

namespace Observable {
// Potential observable
// Takes the real potential stored in state and
// (i) PERIODIC BOUNDARY CONDITIONS:
// =================================
// passes state.V through.
//
// (ii) Homogenous Dirichlet Conditions:
// ====================================
// converts it into a (2*box.N + 2) real vector of the form
//
// [ 1/(4π box.L)            ] <- from boundary condition
// [ V(-box.xmax)            ]
// [ V(-box.xmax + box.dx)   ]
// [        ...              ]
// [ V(-box.dx)              ]
// [ V(0)                    ] <- from regularity (see below)
// [ V(box.dx)               ]
// [        ...              ]
// [ V(box.xmax-box.dx)      ]
// [ V(box.xmax)             ]
//
// i.e. we use the nonredundant half of the radial problem to extend V to a
// symmetric array. This allows other observables to use one common
// implementation for both boundary condition types by reusing the
// Observable::Potential (at the expense of twice as much memory and runtime)

class Potential : public ObservableFunctor {
    const Domain box;
    double t_prev;
    fftw_plan_ptr dst;
    DynamicVector<double> V;

   public:
    Potential(const Parameters& p, const Cosmology&)
        : box(p), t_prev(-1), dst(nullptr), V(0, 0) {
        if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
            V.resize(2 * box.N + 2);
            auto ptr = V.data() + (box.N + 2);
            dst = make_fftw_plan_r2r_1d(box.N, ptr, ptr, FFTW_RODFT00,
                                        FFTW_ESTIMATE);
        }
    }
    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
                // Homogenenoues Dirichlet Conditions at r=0 are artificial in
                // the radial problem and arise due to the transformation
                // V(r) = 4π r V_true(r)
                //
                // Due to the regularity condition of V_true at the origin,
                // i.e.:
                //
                // ∂_r V_true(r) = 0
                //
                // we have:
                //
                // V_true(r) = state.V                          if r>0
                //           = 1/(4π) * ∂_r(4π r state.V)       if r=0
                //
                // which we compute by a spectral derivative.

                auto r = linspace(box.N, box.xmin, box.xmax);
                auto V_r = subvector(V, box.N + 2, box.N);
                V_r = 4 * M_PI * r * state.V;
                fftw_execute(dst.get());
                auto k = linspace(box.N, box.kmin, box.kmax);
                V_r *= k / (box.N + 1);
                V[box.N + 1] = 1.0 / (4 * M_PI) * sum(V_r);
                V_r = state.V;
                auto V_mr = subvector(V, 1, box.N);
                V_mr = reverse(V_r);
                V[0] = 1 / (4 * M_PI * box.L);
            }
        }
        if (box.bc == Domain::BoundaryCondition::Periodic) {
            return state.V;
        }
        return V;
    }

    REGISTER(Potential)
};

}  // namespace Observable
