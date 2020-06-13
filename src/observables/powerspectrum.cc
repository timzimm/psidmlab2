#include <cstddef>
#include "blaze/math/dense/DenseVector.h"
#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"

namespace Observable {

// Computes the single realization matter power spectrum
// P(k) = 1/L * |delta_k|^2
class MatterPowerSpectrum : public ObservableFunctor {
    const Domain box;
    double t_prev;
    DynamicVector<double> P;
    DynamicVector<std::complex<double>> delta_k;
    fftw_plan_ptr r2c;

   public:
    MatterPowerSpectrum(const Parameters& p, const Cosmology&)
        : box(p), t_prev(-1), P(box.N), delta_k(box.N / 2 + 1), r2c(nullptr) {
        auto out = reinterpret_cast<fftw_complex*>(delta_k.data());
        r2c = make_fftw_plan_dft_r2c(box.N, P.data(), out, FFTW_ESTIMATE);
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&)
        override {
        if (t_prev < state.tau) {
            size_t N2 = box.N;
            N2 *= box.N;
            P.resize(box.N);
            P = delta_from(state);
            auto out = reinterpret_cast<fftw_complex*>(delta_k.data());
            fftw_execute_dft_r2c(r2c.get(), P.data(), out);
            P = box.L_phys / N2 * real(conj(delta_k) * delta_k);
            P.resize(box.N / 2 + 1);
        }
        return P;
    }

    REGISTER(MatterPowerSpectrum)
};

}  // namespace Observable
