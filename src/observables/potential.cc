#include "observables_common.h"

namespace Observable {

class Potential : public ObservableFunctor {
   public:
    Potential(const Parameters& p, const Cosmology&) {}
    ReturnType compute(
        const SimState& state,
        const std::unordered_map<std::string,
                                 std::unique_ptr<ObservableFunctor>>& obs) {
        return state.V;
    }

    REGISTER(Potential)
};

}  // namespace Observable
