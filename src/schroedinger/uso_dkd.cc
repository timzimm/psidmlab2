#include "schroedinger/uso_dkd.h"
#include <cassert>
#include "cosmology.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

USO_DKD::USO_DKD(const Parameters& p)
    : cosmo{p},
      pot{PotentialMethod::make(p.get<std::string>("potential"), p)},
      N{p.get<size_t>("N")},
      L{p.get<double>("L")},
      K(N, N),
      D(N, N),
      wavenumbers(N),
      forwards(nullptr),
      backwards(nullptr) {
    for (int i = 0; i < N; ++i) {
        double k = 2 * M_PI / L * i;
        wavenumbers[i] = k * k;
    }
    // TODO: Maybe one value is enough??
    CCV fft_dummy(N);
    auto in = reinterpret_cast<fftw_complex*>(fft_dummy.data());
    forwards = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void USO_DKD::transform_matrix(const fftw_plan& plan, CCM& matrix_in,
                               CCM& matrix_out) {
    assert(matrix_in.columns() == matrix_out.columns());

    for (int i = 0; i < matrix_in.columns(); ++i) {
        fftw_complex* row_in =
            reinterpret_cast<fftw_complex*>(&(matrix_in(0, i)));
        fftw_complex* row_out =
            reinterpret_cast<fftw_complex*>(&(matrix_out(0, i)));
        fftw_execute_dft(plan, row_in, row_out);
    }
}

// Kick Operator - returns a blaze expression
auto USO_DKD::kick(const CCM& psis_in_k, const double dt, const double weight) {
    auto diag_K = blaze::diagonal(K);
    diag_K = blaze::exp(-1.0 * weight / 2 * cmplx(0, 1) * wavenumbers *
                        wavenumbers * dt);

    return K * psis_in_k;
}

// Drift Operator - returns a blaze expression
auto USO_DKD::drift(const CCM& psis_in_x, const RCV& V, const double dt,
                    const double t, const double weight) {
    auto diag_D = blaze::diagonal(D);
    diag_D =
        blaze::exp(-1.0 * weight * cmplx(0, 1) * cosmo.a_of_tau(t) * V * dt);
    // This is a dense-sparse product
    return D * psis_in_x;
}

USO_DKD::~USO_DKD() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
}

void USO_DKD::step(SimState& state) {
    double dt = state.dtau;
    double t = state.tau;
    CCM& psis = state.psis;
    RCV& V = state.V;

    // Drift step
    // It is crucial to always use the most up-to date version of psi.
    // Therefore, we cannot reuse the potential of the last step but
    // have to recalculate it.
    pot->solve(state);
    psis = drift(psis, V, t, dt, 1.0 / 2);

    transform_matrix(forwards, psis, psis);

    psis = kick(psis, dt, 1.0);

    transform_matrix(backwards, psis, psis);

    pot->solve(state);
    psis = drift(psis, V, t, dt, 1.0 / 2);

    // psis and V are now @ tau + dtau
    // Update time information
    state.tau += dt;
    state.a = cosmo.a_of_tau(t + dt);
}
}  // namespace Schroedinger
