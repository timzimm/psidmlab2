#include "evolution/bfsp.h"
#include "H5Opublic.h"
#include "blaze/math/dense/DynamicVector.h"
#include "cosmology.h"
#include "fftw.h"
#include "fftw3.h"
#include "parameters.h"
#include "state.h"

#include <cmath>

namespace DNGF {
using namespace blaze;

BFSP::BFSP(const Parameters& p, const SimState& state, const Cosmology& cosmo_)
    : DefaultDriver(p),
      cosmo(cosmo_),
      box(p),
      phi(box.N),
      bmin(0),
      bmax(0),
      alpha(0),
      // Normalization
      norm0(std::sqrt(sum(box.dx * rho_from(state)))),
      dst(nullptr),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)} {
    if (box.bc != Domain::BoundaryCondition::HomogeneousDirichlet) {
        std::cout << ERRORTAG("Wrong Boundary Conditions") << std::endl;
        exit(1);
    }
    std::cout << INFOTAG("Generating transformation plan for Gradient Descent")
              << std::endl;
    dst = make_fftw_plan_r2r_1d(box.N, phi.data(), phi.data(), FFTW_RODFT00,
                                FFTW_MEASURE);
}

// As stated in W.Bao Journal of Compuatational Physics 219 (2006) 836-854
double BFSP::next_dt(const SimState& state) const {
    const double security_factor = 0.9;

    bmin = min(abs(cosmo.a_of_tau(state.tau) * state.V));
    bmax = max(abs(cosmo.a_of_tau(state.tau) * state.V));
    alpha = 0.5 * (bmin + bmax);
    return security_factor * 2.0 / (bmin + bmax);
}

void BFSP::step(SimState& state, const double dt) const {
    auto& aVphi = state.V;

    state.transform(SimState::Representation::Position);
    // Ground state is real
    phi = real(state.psi);

    // Potential was already updated in last step. See below
    aVphi = (alpha - cosmo.a_of_tau(state.tau) * state.V) * phi;

    fftw_execute_r2r(dst.get(), aVphi.data(), aVphi.data());
    fftw_execute_r2r(dst.get(), phi.data(), phi.data());

    auto k = linspace(box.N, box.kmin, box.kmax);

    // Forward Backward Euler step + optimal stability factor of next_dt()
    phi = 2.0 / (2 + dt * (2 * alpha + k * k)) * (phi + dt * aVphi);

    // Return to real space
    fftw_execute_r2r(dst.get(), phi.data(), phi.data());

    // Renormalization to norm0 (takes also care of unormalized DST)
    phi *= norm0 / std::sqrt(box.dx * sum(phi * phi));

    // Retransfer to complex state and compute potential for next step
    state.psi = phi;
    pot->solve(state);

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
}

}  // namespace DNGF
