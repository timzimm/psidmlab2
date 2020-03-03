#include "observables_common.h"

namespace Observable {

// Wave function observable
// Takes the complex wavefunction stored in state and converts it into a
// N x 2 real matrix
class WaveFunction : public ObservableFunctor {
    double t_prev;
    DynamicMatrix<double> psi_re_im;

   public:
    WaveFunction(const Parameters& p, const Cosmology&)
        : t_prev(-1), psi_re_im{p["Simulation"]["N"], 2} {};

    ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            column(psi_re_im, 0) = real(state.psi);
            column(psi_re_im, 1) = imag(state.psi);
        }
        return psi_re_im;
    }

    REGISTER(WaveFunction)
};

}  // namespace Observable
