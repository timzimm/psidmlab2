#include "observables_common.h"

namespace Observable {

class Entropy : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const int N;
    const double dx;
    const double dk;
    double t_prev;
    DynamicVector<double> S;
    std::unique_ptr<ObservableFunctor> phasespace;

   public:
    Entropy(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_),
          cosmo(cosmo_),
          N{p["Simulation"]["N"]},
          dx{p["Simulation"]["L"].get<double>() / N},
          dk(M_PI / dx),
          t_prev(-1),
          S(1),
          phasespace(nullptr) {}

    ReturnType compute(
        const SimState& state,
        const std::unordered_map<std::string,
                                 std::unique_ptr<ObservableFunctor>>& obs) {
        // rho could be either husimi (row-major) or wigner (column-major).
        auto boltzmann_entropy = [&](const auto& rho) {
            // both x and k are periodic so trapezodial integration is very
            // accurate.
            return -dx * dk * sum(rho * log(rho));
        };
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // rho already known?
            if (phasespace == nullptr) {
                // Check if there exists a PhaseSpaceDistribution observable
                const std::string name = "Observable::PhaseSpaceDistribution";
                auto result = obs.find(name);
                if (result == obs.end()) {
                    phasespace = ObservableFunctor::make(name, p, cosmo);
                } else {
                    auto rho = result->second->compute(state, obs);
                    S[0] = boost::apply_visitor(boltzmann_entropy, rho);
                }
            }
            auto rho = phasespace->compute(state, obs);
            S[0] = boost::apply_visitor(boltzmann_entropy, rho);
        }
        return S;
    }
    REGISTER(Entropy)
};

}  // namespace Observable
