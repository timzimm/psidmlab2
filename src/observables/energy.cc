#include "blaze/math/ReductionFlag.h"
#include "cosmology.h"
#include "observables_common.h"

namespace Observable {

class EnergyDensity : public ObservableFunctor {
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    // Sparse matrix holding the finite difference approximation
    // of the first derivative
    CompressedMatrix<double> grad;
    // 3 x N matrix of energy densities:
    // (i)   kinetic energy
    // (ii)  potential energy
    // (iii) virial density
    DynamicMatrix<double> energies;

   public:
    EnergyDensity(const Parameters& p, const Cosmology& cosmo_)
        : cosmo(cosmo_),
          box{p},
          t_prev(-1),
          grad(box.N, box.N + 4),
          energies(3, box.N) {
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

            const double a = cosmo.a_of_tau(state.tau);

            // (i) Kinetic energy density
            auto E_kinetic = row(energies, 0);
            // Cyclic extension of psi to allow gradient computation across
            // boundary
            auto psi = elements(state.psi, cyclic_extension, box.N + 4);
            // Sparse matrix-vector product
            auto nabla_psi = grad * psi;
            E_kinetic =
                trans(0.5 / (a * a) * real(conj(nabla_psi) * nabla_psi));

            // (ii) Potential energy density
            auto E_pot = row(energies, 1);
            E_pot =
                trans(0.5 / a * state.V * real(conj(state.psi) * state.psi));

            // (iii) Virial
            // Assuming V to be arbitrary (i.e. potentially non homogeneous) we
            // compute the virial in its general form, i.e. x grad V |psi|^2
            auto virial = row(energies, 2);
            auto V = elements(state.V, cyclic_extension, box.N + 4);
            auto x = linspace(box.N, box.xmin, box.xmax);
            virial = trans(1.0 / a * x * (grad * V) *
                           real(conj(state.psi) * state.psi));
        }

        return energies;
    }
    REGISTER(EnergyDensity)
};

class Energy : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    DynamicVector<double> energies;

   public:
    Energy(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_), cosmo(cosmo_), box(p), t_prev(-1), energies(3, 0) {}

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // Check if there exists a EnergyDensity observable
            const std::string name = "EnergyDensity";
            const std::string full_name = "Observable::" + name;
            auto result = obs.find(name);
            if (result == obs.end()) {
                // Allocate new EntropyDensity observable and
                // store it in global obs map
                obs[name] = ObservableFunctor::make(full_name, p, cosmo);
            }
            ObservableFunctor* edensity = obs[name].get();
            ReturnType energy_density = edensity->compute(state, obs);
            // Trapezodial integration in x-space

            energies =
                box.dx * sum<rowwise>(boost::get<const DynamicMatrix<double>&>(
                             energy_density));
        }

        return energies;
    }
    REGISTER(Energy)
};

}  // namespace Observable
