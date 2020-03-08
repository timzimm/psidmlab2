#include "interfaces.h"
#include "observables_common.h"

namespace Observable {

class EntropyDensity : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const bool husimi;
    const int N;
    const double dk;
    double t_prev;
    DynamicVector<double> s;

   public:
    EntropyDensity(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_),
          cosmo(cosmo_),
          husimi{p["Observables"]["PhaseSpaceDistribution"]["sigma_x"] > 0},
          N{p["Simulation"]["N"]},
          dk{2 * M_PI / p["Simulation"]["L"].get<double>()},
          t_prev(-1),
          s(1) {}

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        // rho could either be husimi (row-major) or wigner (column-major).
        // which are different types.
        auto wehrl_entropy_density = [&](const auto& rho) {
            // k is periodic (because x is discrete) so trapezodial integration
            // is very accurate.
            // We collapse the matrix columnwise. The result is a row vector.
            // Transpose it to get a column vector
            return -dk * trans(sum<columnwise>(rho % log(rho)));
        };
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // Check if there exists a PhaseSpaceDistribution observable
            const std::string name = "PhaseSpaceDistribution";
            auto result = obs.find(name);
            if (result == obs.end()) {
                // Allocate new PhaseSpaceDistribution observable and
                // store it in global obs map
                obs[name] = ObservableFunctor::make(name, p, cosmo);
            }
            ObservableFunctor* phasespace = obs[name].get();
            ReturnType rho = phasespace->compute(state, obs);
            // ReturnType is a variant that can also hold refs to vectors.
            // sum<columswise> only operates on matrices, so we have to
            // explicitly extract the type.
            if (husimi) {
                s = wehrl_entropy_density(
                    boost::get<const DynamicMatrix<double>&>(rho));
            } else {
                s = wehrl_entropy_density(
                    boost::get<const DynamicMatrix<double, columnMajor>&>(rho));
            }
        }
        return s;
    }
    REGISTER(EntropyDensity)
};

class Entropy : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const int N;
    const double dx;
    double t_prev;
    DynamicVector<double> S;

   public:
    Entropy(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_),
          cosmo(cosmo_),
          N{p["Simulation"]["N"]},
          dx{p["Simulation"]["L"].get<double>() / N},
          t_prev(-1),
          S(1) {}

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        // Trapezodial integration in x-space
        auto integrate_density = [&](const auto& entropy_density) {
            return dx * sum(entropy_density);
        };
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // Check if there exists a EntropyDensity observable
            const std::string name = "EntropyDensity";
            const std::string full_name = "Observable::" + name;
            auto result = obs.find(name);
            if (result == obs.end()) {
                // Allocate new EntropyDensity observable and
                // store it in global obs map
                obs[name] = ObservableFunctor::make(full_name, p, cosmo);
            }
            ObservableFunctor* sdensity = obs[name].get();
            ReturnType s = sdensity->compute(state, obs);
            S[0] = boost::apply_visitor(integrate_density, s);
        }
        return S;
    }
    REGISTER(Entropy)
};

}  // namespace Observable
