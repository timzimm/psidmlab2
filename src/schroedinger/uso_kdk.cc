#include "schroedinger/uso_kdk.h"
#include <cassert>
#include "cosmology.h"
#include "parameters.h"
#include "state.h"

namespace Schroedinger {

USO_KDK::USO_KDK(const Parameters& p, const Cosmology& cosmo_)
    : cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      N{p["Simulation"]["N"].get<size_t>()},
      L{p["Simulation"]["L"].get<double>()},
      firstStep(true),
      K(N, N),
      D(N, N),
      wavenumbers(N),
      forwards(nullptr),
      backwards(nullptr),
      forwards_op(nullptr),
      backwards_op(nullptr),
      psis_cached(N, N) {
    for (int i = 0; i < N; ++i) {
        double k = 2 * M_PI / L * i;
        wavenumbers[i] = k * k;
    }
    CCV fft_dummy(N);
    auto in = reinterpret_cast<fftw_complex*>(fft_dummy.data());
    auto out = reinterpret_cast<fftw_complex*>(psis_cached.data());
    forwards = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    forwards_op = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    backwards_op = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
}

void USO_KDK::transform_matrix(const fftw_plan& plan, CCM& matrix_in,
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
auto USO_KDK::kick(const CCM& psis_in_k, const double dt, const double weight) {
    auto diag_K = blaze::diagonal(K);
    diag_K = blaze::exp(-1.0 * weight / 2 * cmplx(0, 1) * wavenumbers *
                        wavenumbers * dt);

    // psi_i is the i-the row in the psis matrix. Therfore the
    // multiplictaion order changes.
    // This is a dense-sparse product
    return K * psis_in_k;
}

// Drift Operator - returns a blaze expression
auto USO_KDK::drift(const CCM& psis_in_x, const RCV& V, const double dt,
                    const double t, const double weight) {
    auto diag_D = blaze::diagonal(D);
    diag_D =
        blaze::exp(-1.0 * weight * cmplx(0, 1) * cosmo.a_of_tau(t) * V * dt);
    // This is a dense-sparse product
    return D * psis_in_x;
}

USO_KDK::~USO_KDK() {
    fftw_destroy_plan(forwards);
    fftw_destroy_plan(backwards);
    fftw_destroy_plan(forwards_op);
    fftw_destroy_plan(backwards_op);
}

void USO_KDK::step(SimState& state) {
    double dt = state.dtau;
    double t = state.tau;
    CCM& psis = state.psis;
    RCV& V = state.V;

    // We can spare the initial FFT if we use the cached psi_in_k
    // representation of the last step
    if (firstStep) {
        transform_matrix(forwards_op, psis, psis_cached);
        firstStep = false;
    }

    psis = kick(psis_cached, dt, 1.0 / 2);

    transform_matrix(backwards, psis, psis);

    // Recompute potential
    pot->solve(state);
    psis = drift(psis, V, t, dt, 1.0);

    transform_matrix(forwards, psis, psis);

    // Update cached psis
    psis_cached = kick(psis, dt, 1.0 / 2);

    transform_matrix(backwards_op, psis_cached, state.psis);

    // psis and V are now @ tau + dtau
    // Update time information
    state.tau += dt;
    state.a = cosmo.a_of_tau(t + dt);
}

}  // namespace Schroedinger
