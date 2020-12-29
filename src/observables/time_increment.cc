#include "observables_common.h"

namespace Observable {

// Time Increment observable
class TimeIncrement : public ObservableFunctor {
    double t_prev;
    DynamicVector<double> dt;

  public:
    TimeIncrement(const Parameters &p, const Cosmology &) : t_prev(-1), dt(1){};

    ReturnType
    compute(const SimState &state,
            std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>
                &obs) override {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            dt[0] = state.dtau;
        }
        return dt;
    }

    REGISTER(TimeIncrement)
};

} // namespace Observable
