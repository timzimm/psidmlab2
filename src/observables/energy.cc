#include "cosmology.h"
#include "observables_common.h"

namespace Observable {

class Energy : public ObservableFunctor {
    const Cosmology& cosmo;
    const int N;
    const double L;
    const double dx;
    // Sparse matrix holding the finite difference approximation
    // of the first derivative
    CompressedMatrix<double> grad;
    DynamicVector<double> energies;

   public:
    Energy(const Parameters& p, const Cosmology& cosmo_)
        : cosmo(cosmo_),
          N{p["Simulation"]["N"]},
          L{p["Simulation"]["L"]},
          dx(L / N),
          grad(N, N + 4),
          energies(4, 0) {
        grad.reserve(4 * N);

        // 5 point stencil for first derivative
        for (int i = 0; i < N; ++i) {
            grad.append(i, i, 1.0 / (12 * dx));
            grad.append(i, i + 1, -2.0 / (3 * dx));
            grad.append(i, i + 3, 2.0 / (3 * dx));
            grad.append(i, i + 4, -1.0 / (12 * dx));
            grad.finalize(i);
        }
    };

    ReturnType compute(
        const SimState& state,
        const std::unordered_map<std::string,
                                 std::unique_ptr<ObservableFunctor>>& obs) {
        using namespace blaze;
        auto cyclic_extension = [N = N](int i) { return (N - 2 + i) % N; };

        const double a = cosmo.a_of_tau(state.tau);

        // Kinetic Energy via trapezodial rule
        double& E_kinetic = energies[0];
        auto psi = elements(state.psi, cyclic_extension, N + 4);
        auto nabla_psi = grad * psi;
        E_kinetic = 0.5 * dx / (a * a) * sum(real(conj(nabla_psi) * nabla_psi));

        // Potential Energy via trapezodial rule
        double& E_pot = energies[1];
        E_pot =
            0.5 * dx / a * sum(state.V * (real(conj(state.psi) * state.psi)));

        // Total Energy
        double& E_tot = energies[2];
        E_tot = E_kinetic + E_pot;

        // Assuming V to be arbitrary (i.e. potentially non homogeneous) we
        // compute the virial in its general form, i.e. <x grad V>
        double& x_grad_V = energies[3];
        auto V = elements(state.V, cyclic_extension, N + 4);
        auto nabla_V = grad * V;
        auto x = linspace(N, -L / 2, L / 2 - dx);
        x_grad_V =
            dx / a * sum(x * nabla_V * (real(conj(state.psi) * state.psi)));

        return energies;
    }
    REGISTER(Energy)
};

}  // namespace Observable
