#include "interaction/poisson_fd.h"
#include "blaze_utils.h"
#include "parameters.h"

namespace Poisson {
FD::FD(const Parameters& p)
    : dx{p["Simulation"]["L"].get<double>() / p["Simulation"]["N"].get<int>()},
      N{p["Simulation"]["N"].get<size_t>()},
      dl(N, 1),
      d(N, -2),
      du(N, 1),
      du2(N, 0),
      ipiv(N) {
    // The LU decomposition never changes for Poissons equation. Therefore
    // we compute it once and reuse it every time we need the potential for
    // a new source term.
    gctrf(dl, d, du, du2, ipiv);
}

void FD::solve(SimState& state) { solve(state.V, delta_from(state)); }

void FD::solve(blaze::DynamicVector<double>& V,
               const blaze::DynamicVector<double>& source) {
    using namespace blaze;
    // Calculate source term as matrix
    RCM source_mat = expand(dx * dx * source, 1);

    gctrs(dl, d, du, du2, ipiv, source_mat);

    // Solution is only determined up to an additive constant. We fix it by
    // requiring a vanishing mean. This is in accordance with the behavior
    // of Poisson::FFT, where the DC mode is set to 0 explicitly
    double mean = 1.0 / N * sum(source_mat);
    V = column(map(source_mat, [&mean](double V) { return V - mean; }), 0);
}

}  // namespace Poisson
