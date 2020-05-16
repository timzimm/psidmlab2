#include <limits>
#include "blaze/math/ReductionFlag.h"
#include "blaze/system/TransposeFlag.h"
#include "cosmology.h"
#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"
#include "state.h"

using namespace std::complex_literals;

namespace Observable {

// Computation of the forces:
// (i)   F_G = - ∂_x V                                        (gravity)
// (ii)  F_Q = 1/2 * ∂_x [ ∂^2_x{ sqrt(rho) } / sqrt(rho) ]   (quantum force)

class Force : public ObservableFunctor {
    const Cosmology& cosmo;
    const Domain box;
    double t_prev;
    // Sparse matrix holding the finite difference approximation
    // of the first derivative
    CompressedMatrix<double> grad;

    // 2 x N matrix of forces:
    DynamicMatrix<double> forces;
    DynamicVector<std::complex<double>> gradient;
    fftw_plan_ptr c2r, r2c;

   public:
    Force(const Parameters& p, const Cosmology& cosmo_)
        : cosmo(cosmo_),
          box{p},
          t_prev(-1),
          grad(box.N, box.N + 4),
          gradient(box.N / 2 + 1),
          forces(2, box.N) {
        auto cmplx_p = reinterpret_cast<fftw_complex*>(gradient.data());
        auto real_p = forces.data();
        r2c = make_fftw_plan_dft_r2c(box.N, real_p, cmplx_p, FFTW_ESTIMATE);
        c2r = make_fftw_plan_dft_c2r(box.N, cmplx_p, real_p, FFTW_ESTIMATE);
        grad.reserve(4 * box.N);
        // 5 point stencil for first derivative
        for (int i = 0; i < box.N; ++i) {
            grad.append(i, i, 1.0 / (12 * box.dx));
            grad.append(i, i + 1, -2.0 / (3 * box.dx));
            grad.append(i, i + 3, 2.0 / (3 * box.dx));
            grad.append(i, i + 4, -1.0 / (12 * box.dx));
            grad.finalize(i);
        }
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (t_prev < state.tau) {
            t_prev = state.tau;

            // k-grid
            auto k0p = linspace(box.N / 2 + 1, 0.0, box.N / 2 * box.dk);

            // aliases in k-space
            auto gradient_plus = subvector(gradient, 0, box.N / 2 + 1);

            auto gradient_p = reinterpret_cast<fftw_complex*>(gradient.data());

            // (ii) Quantum pressure force
            auto sqrt_rho = row(forces, 0);
            auto F_Q_wip = row(forces, 1);
            sqrt_rho = trans(sqrt(rho_from(state)));
            fftw_execute_dft_r2c(r2c.get(), sqrt_rho.data(), gradient_p);
            gradient_plus *= -k0p * k0p;
            fftw_execute_dft_c2r(c2r.get(), gradient_p, F_Q_wip.data());
            F_Q_wip /= 2 * box.N * sqrt_rho;

            // Since sqrt_rho is small F_Q_wip explodes and becomes very spiky.
            // Taking spectral derivatives then becomes unfavourable due to
            // ringing. That's whay the last derivative is computed via a finite
            // difference approximation

            // Cyclic extension of F_Q to allow gradient computation across
            // boundary
            auto cyclic_extension = [N = box.N](int i) {
                return (N - 2 + i) % N;
            };
            auto F_Q_pi = elements(F_Q_wip, cyclic_extension, box.N + 4);
            // Sparse product
            row(forces, 0) = F_Q_pi * trans(grad);

            // (ii) gravitational force
            auto F_G = row(forces, 1);
            F_G = trans(state.V);
            fftw_execute_dft_r2c(r2c.get(), F_G.data(), gradient_p);
            gradient_plus *= 1.0i * k0p;
            fftw_execute_dft_c2r(c2r.get(), gradient_p, F_G.data());
            F_G /= -1.0 * box.N;
        }

        return forces;
    }
    REGISTER(Force)
};

}  // namespace Observable
