#include "observables_common.h"

namespace Observable {

// Wave function observable
// Takes the complex wavefunction stored in state and converts it into a
// N x 2 real matrix
class WaveFunction : public ObservableFunctor {
    const Domain box;
    double t_prev;
    DynamicMatrix<double> psi_re_im;

   public:
    WaveFunction(const Parameters& p, const Cosmology&)
        : box(p), t_prev(-1), psi_re_im(box.N, 2){};

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) override {
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
