#include <limits>
#include "cosmology.h"
#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"
#include "state.h"

using namespace std::complex_literals;

namespace Observable {

// Computation of the peculiar velocity:

class PeculiarVelocity : public ObservableFunctor {
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;

    DynamicVector<double> v;
    DynamicVector<std::complex<double>> gradient;
    fftw_plan_ptr forward, backward;

    // Convert code velocity to km/s
    void convert_to_physical_units(const SimState& state) {
        static const double c = 2.997 * 1e8;                    // m/s
        static const double hbar = 6.636 / (2 * M_PI) * 1e-34;  // Js
        static const double m_per_Mpc = 3.2407 * 1e-23;
        static const double eV_per_J = 1.6021 * 1e-19;
        static const double L = box.L_phys / box.L;  // Mpc/code length
        static const double hbar_over_m = (c * c * hbar) /
                                          (box.m22 * 1e-22 * eV_per_J) * 1e-3 *
                                          m_per_Mpc;  // Mpc km /s
        v *= hbar_over_m / L;
    }

   public:
    PeculiarVelocity(const Parameters& p, const Cosmology& cosmo_)
        : cosmo(cosmo_), box{p}, t_prev(-1), gradient(box.N), v(box.N) {
        auto ptr = reinterpret_cast<fftw_complex*>(gradient.data());
        forward =
            make_fftw_plan_dft(box.N, ptr, ptr, FFTW_FORWARD, FFTW_ESTIMATE);
        backward =
            make_fftw_plan_dft(box.N, ptr, ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // positive and negative k-multipliers combined with FFT
            // normalization
            auto k0p = linspace(box.N / 2 + 1, 0.0,
                                1.0 / box.N * (box.N / 2 * box.dk));
            const int NN = (box.N % 2) ? box.N + 1 : box.N;
            auto km = linspace(NN / 2 - 1, 1.0 / box.N * (-NN / 2 + 1) * box.dk,
                               1.0 / box.N * (-1 * box.dk));

            // aliases in k-space
            auto gradient_plus = subvector(gradient, 0, box.N / 2 + 1);
            auto gradient_minus =
                subvector(gradient, box.N / 2 + 1, NN / 2 - 1);

            gradient = state.psi;
            auto ptr = reinterpret_cast<fftw_complex*>(gradient.data());
            fftw_execute_dft(forward.get(), ptr, ptr);
            gradient_plus *= 1.0i * k0p;
            gradient_minus *= 1.0i * km;
            fftw_execute_dft(backward.get(), ptr, ptr);
            v = real(gradient);

            const double a = cosmo.a_of_tau(state.tau);
            v = 1.0 / (a * rho_from(state)) * imag(gradient * conj(state.psi));
            if (box.physical) convert_to_physical_units(state);
        }

        return v;
    }
    REGISTER(PeculiarVelocity)
};

}  // namespace Observable
