#include "blaze/math/ReductionFlag.h"
#include "cosmology.h"
#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"
#include "state.h"

using namespace std::complex_literals;

namespace Observable {

// For a NLSE of the form
//
// i∂_t psi = -1/2∂_xx psi + a (G * |psi|^2) psi        (*)
//
// and under PERIODIC BOUNDARY CONDITIONS the generating energy functional is:
//          xmax                      xmax
// H[psi] = ∫dx[ 1/2 |∂_x psi|^2 + a/2∫ dx' G(x,x') |psi(x')|^2 |psi(x)|^2 ]
//          xmin                      xmin
//        = <T> + a/2 <V>
//
// If a=const, H is a conserved quantity of (*). Note the 1/2 to account for the
// SYMMETRIC kernel G(x,x') = G(x',x)
//
// The class EnergyDensity computes the NON-INTEGRATED (energy densities)
// quantities:
//
//  energies[0] = T         = 1/2 |∂_x psi|^2
//  energies[1] = V         = ∫ dx' G(x,x') |psi(x')|^2 |psi(x)|^2
//  energies[2] = virial    = x ∂_x( ∫ dx' G(x,x') |psi(x')|^2 ) |psi|^2
//
//  by a finite difference approximation.

class EnergyDensity : public ObservableFunctor {
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    // 3 x N matrix of energy densities:
    // (i)   kinetic energy
    // (ii)  potential energy
    // (iii) virial density
    DynamicMatrix<double> energies;
    DynamicVector<std::complex<double>> gradient;
    fftw_plan_ptr c2c_f, c2c_b, c2r, r2c;

   public:
    EnergyDensity(const Parameters& p, const Cosmology& cosmo_)
        : cosmo(cosmo_),
          box{p},
          t_prev(-1),
          gradient(box.N),
          energies(3, box.N) {
        auto cmplx_p = reinterpret_cast<fftw_complex*>(gradient.data());
        auto real_p = &energies(2, 0);
        c2c_f = make_fftw_plan_dft(box.N, cmplx_p, cmplx_p, FFTW_FORWARD,
                                   FFTW_ESTIMATE);
        c2c_b = make_fftw_plan_dft(box.N, cmplx_p, cmplx_p, FFTW_BACKWARD,
                                   FFTW_ESTIMATE);
        r2c = make_fftw_plan_dft_r2c(box.N, real_p, cmplx_p, FFTW_ESTIMATE);
        c2r = make_fftw_plan_dft_c2r(box.N, cmplx_p, real_p, FFTW_ESTIMATE);
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

            // (i) Kinetic energy density
            auto E_kinetic = row(energies, 0);
            gradient = state.psi;
            fftw_execute(c2c_f.get());
            gradient_plus *= 1.0i * k0p;
            gradient_minus *= 1.0i * km;
            fftw_execute(c2c_b.get());
            E_kinetic = trans(0.5 * real(conj(gradient) * gradient));

            // (ii) Potential energy density
            auto E_pot = row(energies, 1);
            E_pot = trans(state.V * rho_from(state));

            // (iii) Virial
            auto virial = row(energies, 2);
            virial = trans(state.V);
            fftw_execute(r2c.get());
            gradient_plus *= 1.0i * k0p;
            fftw_execute(c2r.get());

            auto x = linspace(box.N, box.xmin, box.xmax);
            virial *= trans(x * rho_from(state));
        }

        return energies;
    }
    REGISTER(EnergyDensity)
};

// The class Energy integrates the quanitites computed by EnergyDensity
//
//  energies[0] = <T>       = ∫dx 1/2 |∂_x psi|^2
//  energies[1] = <V>       = ∫dx ∫dx' G(x,x') |psi(x')|^2 |psi(x)|^2
//  energies[2] = <virial>  = ∫dx x ∂_x( ∫ dx' G(x,x') |psi(x')|^2 ) |psi|^2

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
                // Allocate new EnergyDensity observable and
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
