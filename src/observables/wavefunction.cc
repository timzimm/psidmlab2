#include "domain.h"
#include "fftw.h"
#include "fftw3.h"
#include "observables_common.h"

namespace Observable {

// Wave function observable
// Takes the complex wavefunction stored in state and converts it into:
// (i) PERIODIC BOUNDARY CONDITIONS:
// =================================
// a box.N x 2 real matrix of the form
//
// [ Re( psi(box.xmin) )            Im( psi(box.xmin) )             ]
// [ Re( psi(box.xmin + box.dx) )   Im( psi(box.xmin + box.dx) )    ]
// [            ...                             ...                 ]
// [ Re( psi(box.xmax - box.dx) )   Im( psi(box.xmax - box.dx) )    ]
// [ Re( psi(box.xmax) )            Im( psi(box.xmax) )             ]
//
// (ii) Homogenous Dirichlet Conditions:
// ====================================
// a (box.N + 1) x 2 real matrix of the form
//
// [ Re( psi(0) )                   Im( psi(0) )                    ]
// [ Re( psi(box.dx) )              Im( psi(box.dx) )               ]
// [            ...                             ...                 ]
// [ Re( psi(box.xmax-box.dx) )     Im( psi(box.xmax-box.dx) )      ]
// [ Re( psi(box.xmax) )            Im( psi(box.xmax) )             ]
class WaveFunction : public ObservableFunctor {
    const Domain box;
    double t_prev;
    fftw_plan_ptr dst;
    DynamicMatrix<double> psi_re_im;

   public:
    WaveFunction(const Parameters& p, const Cosmology&)
        : box(p), t_prev(-1), dst(nullptr), psi_re_im(box.N, 2, 0) {
        if (box.bc == Domain::BoundaryCondition::HomogeneousDirichlet) {
            psi_re_im.resize(box.N + 1, 2);
            const int howmany = 2;
            const int n[] = {box.N, box.N};
            const int rank = 1;
            const int idist = 1;
            const int odist = idist;
            const int istride = psi_re_im.spacing();
            const int ostride = istride;
            const fftw_r2r_kind kinds[] = {FFTW_RODFT00, FFTW_RODFT00};
            auto ptr = psi_re_im.data() + istride * 1;
            dst = make_fftw_plan_many_r2r(rank, n, howmany, ptr, nullptr,
                                          istride, idist, ptr, nullptr, ostride,
                                          odist, kinds, FFTW_ESTIMATE);
        }
    };

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) override {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            // Under periodic condition nothing (apart from splitting into
            // real and imaginary part) is left to be done
            if (box.bc == Domain::BoundaryCondition::Periodic) {
                column(psi_re_im, 0) = real(state.psi);
                column(psi_re_im, 1) = imag(state.psi);
            } else if (box.bc ==
                       Domain::BoundaryCondition::HomogeneousDirichlet) {
                auto psi_r_real = subvector(column(psi_re_im, 0), 1, box.N);
                auto psi_r_imag = subvector(column(psi_re_im, 1), 1, box.N);
                psi_r_real = 1.0 / (box.N + 1) * real(state.psi);
                psi_r_imag = 1.0 / (box.N + 1) * imag(state.psi);

                // Homogenenoues Dirichlet Conditions at r=0 are artificial in
                // the radial problem and arise due to the transformation
                // psi(r) = 2 sqrt(π) r psi_true(r)
                //
                // Due to the regularity condition of psi_true at the origin,
                // i.e.:
                //
                // ∂_r psi_true(r) = 0
                //
                // we have:
                //
                // psi_true(r) = 1/(2 sqrt(π)) * psi(r)/r       if r>0
                //             = 1/(2 sqrt(π)) * ∂_r psi(r)     if r=0
                //
                // which we compute by a spectral derivative.
                fftw_execute(dst.get());

                auto r = linspace(box.N, box.xmin, box.xmax);
                auto k = linspace(box.N, box.kmin, box.kmax);
                psi_r_real *= k;
                psi_r_imag *= k;

                auto psi_r = submatrix(psi_re_im, 1, 0, box.N, 2);
                // [ Re( psi(0) )   Im( psi(0) ) ]
                row(psi_re_im, 0) =
                    1.0 / (2 * std::sqrt(M_PI)) * sum<columnwise>(psi_r);
                psi_r_real = 1.0 / (2 * std::sqrt(M_PI) * r) * real(state.psi);
                psi_r_imag = 1.0 / (2 * std::sqrt(M_PI) * r) * imag(state.psi);
            }
        }
        return psi_re_im;
    }

    REGISTER(WaveFunction)
};  // namespace Observable

}  // namespace Observable
