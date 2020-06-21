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

class EnergyDensity : public ObservableFunctor {
    const Parameters& p;
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    // 3 x N matrix of energy densities:
    // (i)   kinetic energy
    // (ii)  potential energy
    // (iii) virial density
    DynamicMatrix<double> energies;
    DynamicVector<std::complex<double>> gradient;
    fftw_plan_ptr c2c_f, c2c_b, c2r, r2c, dst_cmplx, dct_cmplx, dst;

    void compute_periodic_energies(const SimState& state) {
        // positive and negative k-multipliers combined with FFT
        // normalization
        auto k0p =
            linspace(box.N / 2 + 1, 0.0, 1.0 / box.N * (box.N / 2 * box.dk));
        const int NN = (box.N % 2) ? box.N + 1 : box.N;
        auto km = linspace(NN / 2 - 1, 1.0 / box.N * (-NN / 2 + 1) * box.dk,
                           1.0 / box.N * (-1 * box.dk));

        // aliases in k-space
        auto gradient_plus = subvector(gradient, 0, box.N / 2 + 1);
        auto gradient_minus = subvector(gradient, box.N / 2 + 1, NN / 2 - 1);

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
    void compute_radial_energies(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        // k-multipliers
        auto k = linspace(box.N, box.kmin, box.kmax);

        // (i) Kinetic energy density
        auto E_kinetic = row(energies, 0);
        auto dpsi_r = subvector(gradient, 1, box.N);
        std::complex<double>& dpsi_0 = gradient[0];
        std::complex<double>& dpsi_Np1 = gradient[box.N + 1];
        dpsi_r = state.psi / (box.N + 1);
        std::cout << subvector(gradient, 1, 8) << std::endl;
        fftw_execute(dst_cmplx.get());
        std::cout << subvector(gradient, 1, 8) << std::endl;
        dpsi_r *= k;
        /* dpsi_0 = sum(dpsi_r); */
        dpsi_0 = 0;
        dpsi_Np1 = 0;
        std::cout << subvector(gradient, 1, 8) << std::endl;
        fftw_execute(dct_cmplx.get());
        std::cout << subvector(gradient, 1, 8) << std::endl;
        auto grad = subvector(gradient, 0, box.N + 1);
        E_kinetic = trans(0.5 * real(conj(grad) * grad));

        // (ii) Potential energy density
        auto E_pot = row(energies, 1);
        // To compute (ii) we also require the full V vector (including r=0)
        // Luckily psi_0 = 0 and we get away by just knowing V_r = state.V
        auto E_pot_r = subvector(E_pot, 1, box.N);
        double& E_pot_0 = E_pot[0];
        E_pot_0 = 0;
        E_pot_r = trans(state.V * rho_from(state));

        // (iii) Virial
        auto r = linspace(box.N, box.xmin, box.xmax);
        auto virial = row(energies, 2);
        auto virial_r = subvector(virial, 1, box.N);
        auto& virial_0 = virial[0];
        virial_0 = 0;
        virial_r = trans(4 * M_PI * r * state.V / (box.N + 1));
        fftw_execute(dst.get());
        virial_r *= trans(k);
        fftw_execute(dst.get());

        virial_r *= trans(r * rho_from(state));
    };

   public:
    EnergyDensity(const Parameters& p_, const Cosmology& cosmo_)
        : p(p_),
          cosmo(cosmo_),
          box{p},
          t_prev(-1),
          energies(3, box.N),
          gradient(box.N),
          c2c_f(nullptr),
          c2c_b(nullptr),
          c2r(nullptr),
          r2c(nullptr),
          dst_cmplx(nullptr),
          dct_cmplx(nullptr),
          dst(nullptr) {
        if (box.bc == Domain::BoundaryCondition::Periodic) {
            auto grad_p = reinterpret_cast<fftw_complex*>(gradient.data());
            auto energy_p = &energies(2, 0);
            c2c_f = make_fftw_plan_dft(box.N, grad_p, grad_p, FFTW_FORWARD,
                                       FFTW_ESTIMATE);
            c2c_b = make_fftw_plan_dft(box.N, grad_p, grad_p, FFTW_BACKWARD,
                                       FFTW_ESTIMATE);
            r2c =
                make_fftw_plan_dft_r2c(box.N, energy_p, grad_p, FFTW_ESTIMATE);
            c2r =
                make_fftw_plan_dft_c2r(box.N, grad_p, energy_p, FFTW_ESTIMATE);
        } else if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
            energies.resize(3, box.N + 1);
            gradient.resize(box.N + 2);
            const int howmany = 2;
            const int n[] = {box.N, box.N};
            const int n_dct[] = {box.N + 2, box.N + 2};
            const int rank = 1;
            const int idist = 1;
            const int odist = idist;
            const int istride = 2;
            const int ostride = istride;
            const fftw_r2r_kind dst_kinds[] = {FFTW_RODFT00, FFTW_RODFT00};
            const fftw_r2r_kind dct_kinds[] = {FFTW_REDFT00, FFTW_REDFT00};
            auto grad_cmplx_p =
                reinterpret_cast<fftw_complex*>(gradient.data());
            auto grad_real_p =
                reinterpret_cast<double*>(grad_cmplx_p) + istride * 1;
            auto energy_p = &energies(2, 1);
            dst_cmplx = make_fftw_plan_many_r2r(
                rank, n, howmany, grad_real_p, nullptr, istride, idist,
                grad_real_p, nullptr, ostride, odist, dst_kinds, FFTW_ESTIMATE);
            dct_cmplx = make_fftw_plan_many_r2r(
                rank, n_dct, howmany, grad_real_p, nullptr, istride, idist,
                grad_real_p, nullptr, ostride, odist, dct_kinds, FFTW_ESTIMATE);
            dst = make_fftw_plan_r2r_1d(box.N, energy_p, energy_p, FFTW_RODFT00,
                                        FFTW_ESTIMATE);
        }
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            if (box.bc == Domain::BoundaryCondition::Periodic) {
                compute_periodic_energies(state);
            }
            if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
                compute_radial_energies(state, obs);
            }
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
