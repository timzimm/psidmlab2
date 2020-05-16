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
      // box.dx is not applicable here! Usually N is the number of non-redundant
      // elements on a periodic grid. Here N is the number of elements not being
      // a Dirichlet point. Hence the grid spacing is L/(N+1)
      norm0(std::sqrt(sum(box.L / (box.N + 1) * rho_from(state)))),
      dst(nullptr),
      pot{Interaction::make(p["Simulation"]["interaction"]["name"], p, state)} {
    // FFTW only stores non-redundant modes. Using
    //              phi_i  = phi_N+i    (periodicity)
    //              phi_i  = -phi_(N-i) (oddness)
    // This implies only the first n integers with N=2(n+1) are non-redundant.
    // n starts counting from i=1 because i=0 must be zero due to both
    // conditions above. We use the DST to enforce homogeneous Dirichlet at
    // i=0 and i=N+1. To satisfy the second root we must have N, i.e. number
    // of non-redundant elements on the grid, as odd:
    //
    // phi0     phi1    phi2    phi3    phi4    phi5    phi6    phi7
    // phi0       |               |    -phi4   -phi3   -phi2   -phi1
    //            +------ N ------+      ^
    //                                   0
    if (box.N % 2 == 0) {
        std::cerr
            << ERRORTAG(
                   "N must be odd to enforce homogenous Dirchlet condition")
            << std::endl;
        exit(1);
    }
    // Also note that RODFT00 is its own inverse. Hence one plan is enough.
    dst = make_fftw_plan_r2r_1d(box.N, phi.data(), phi.data(), FFTW_RODFT00,
                                FFTW_ESTIMATE);
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
    const int M = box.N + 1;
    auto& aVphi = state.V;

    state.transform(SimState::Representation::Position);
    // Ground state is real
    phi = real(state.psi);

    // Potential was already updated in last step. See below
    aVphi = (alpha - cosmo.a_of_tau(state.tau) * state.V) * phi;

    fftw_execute_r2r(dst.get(), aVphi.data(), aVphi.data());
    fftw_execute_r2r(dst.get(), phi.data(), phi.data());

    // sqrt of second order spectral derivative factor
    auto mu = M_PI / box.L * linspace(box.N, 1, box.N);

    // Forward Backward Euler step + optimal stability factor of next_dt()
    phi = 2.0 / (2 + dt * (2 * alpha + mu * mu)) * (phi + dt * aVphi);

    // Return to real space
    fftw_execute_r2r(dst.get(), phi.data(), phi.data());

    // Renormalization to norm0 (takes also care of unormalized DST)
    phi *= norm0 / std::sqrt(box.L / M * sum(phi * phi));

    // Retransfer to complex state and compute potential for next step
    state.psi = phi;
    pot->solve(state);

    // State is now @ state.tau + dtau
    state.tau += dt;
    state.n += 1;
}

}  // namespace DNGF
