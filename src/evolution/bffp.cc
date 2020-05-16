#include "evolution/bffp.h"
#include "cosmology.h"
#include "fftw.h"
#include "fftw3.h"
#include "parameters.h"
#include "state.h"

#include <cmath>

namespace DNGF {
using namespace blaze;

BFFP::BFFP(const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : DefaultDriver(p),
      box(p),
      cosmo(cosmo_),
      aVphi_k(box.N / 2 + 1),
      bmin(0),
      bmax(0),
      alpha(0),
      norm0(std::sqrt(sum(box.dx * rho_from(state)))),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)},
      r2c(nullptr),
      c2r(nullptr) {
    auto out = reinterpret_cast<fftw_complex*>(aVphi_k.data());
    r2c = make_fftw_plan_dft_r2c(box.N, state.V.data(), out, FFTW_ESTIMATE);
    c2r = make_fftw_plan_dft_c2r(box.N, out, state.V.data(), FFTW_ESTIMATE);
}

// As stated in W.Bao Journal of Compuatational Physics 219 (2006) 836-854
double BFFP::next_dt(const SimState& state) const {
    const double security_factor = 0.9;
    const double a = cosmo.a_of_tau(state.tau);

    bmin = min(abs(a * state.V));
    bmax = max(abs(a * state.V));
    alpha = 0.5 * (bmin + bmax);
    return security_factor * 2.0 / (bmin + bmax);
}

void BFFP::step(SimState& state, const double dt) const {
    state.transform(SimState::Representation::Position);

    // The ground state is real
    auto phi = real(state.psi);
    auto phi_k = subvector(state.psi, 0, box.N / 2 + 1);

    const double a = cosmo.a_of_tau(state.tau);
    // Potential was already updated in last step. See below
    // Pseudo-spectral approach uses multiplication in x-space
    state.V = (alpha - a * state.V) * phi;
    // Switch to k-space. Result is now in aVphi_k
    fftw_execute(r2c.get());

    // Reuse the real array normally holding state.V to store the real ground
    // state phi
    state.V = phi;
    // Switch to k-space. Result is now in phi_k, i.e. the lower half of
    // state.psi
    auto out = reinterpret_cast<fftw_complex*>(phi_k.data());
    fftw_execute_dft_r2c(r2c.get(), state.V.data(), out);

    // sqrt of second order spectral derivative factor
    auto k = 2 * M_PI / box.L * linspace(box.N / 2 + 1, 0, box.N / 2);

    // Forward Backward Euler step + optimal stability factor of next_dt()
    aVphi_k = 2.0 / (2 + dt * (2 * alpha + k * k)) * (phi_k + dt * aVphi_k);

    // Return to real space the result, phi, is now in state.V
    fftw_execute(c2r.get());

    // Renormalization to norm0
    double norm_correction = norm0 / std::sqrt(box.dx * sum(state.V * state.V));
    state.psi = norm_correction * state.V;

    // Compute potential for next step
    pot->solve(state);

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
}

}  // namespace DNGF
