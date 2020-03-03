#include "observables_common.h"

namespace Observable {

class DensityContrast : public ObservableFunctor {
    const double sigma_x;         // spatial smoothing scale
    const bool husimi;            // compute husimi?
    const bool linear;            // do linear convolution?
    const int N;                  // number of spatial gridpoints
    const double dx;              // spatial resolution
    const int N_kernel;           // symmetric 5-sigma_x interval in points
    double t_prev;                // timestamp of last cached observable
    convolution_ws<double> ws;    // holds husimi delta
    DynamicVector<double> delta;  // holds wigner delta

   public:
    DensityContrast(const Parameters& p, const Cosmology&)
        : sigma_x{p["Observables"]["DensityContrast"]["sigma_x"]},
          husimi(sigma_x > 0),
          linear{p["Observables"]["DensityContrast"]["linear_convolution"]},
          N{p["Simulation"]["N"]},
          dx{p["Simulation"]["L"].get<double>() / N},
          N_kernel(2 * floor(5 * sigma_x / dx) + 1),
          t_prev(-1),
          ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
          delta(husimi ? 0 : N) {
        if (husimi) {
            auto& gaussian = ws.kernel_padded;
            std::iota(std::begin(gaussian), std::end(gaussian), -N_kernel / 2);
            gaussian =
                exp(-0.5 * dx * dx * gaussian * gaussian / (sigma_x * sigma_x));
            gaussian /= sum(gaussian);
        }
    }

    ReturnType compute(
        const SimState& state,
        const std::unordered_map<std::string,
                                 std::unique_ptr<ObservableFunctor>>& obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            if (husimi) {
                // Store delta inside convolution workspace
                ws.signal_padded = delta_from(state);
                discrete_convolution(ws);
                return ws.signal_padded;
            }
            // Store delta in dedicated array
            delta = delta_from(state);
        }

        return delta;
    }
    REGISTER(DensityContrast)
};

}  // namespace Observable
