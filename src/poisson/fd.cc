#include "poisson/fd.h"
#include "blaze/math/UniformVector.h"
#include "lapacke_blaze_wrapper.h"
#include "parameters.h"
#include "state.h"

namespace Poisson {
FD::FD(const Parameters& p)
    : dx{p.get<double>("L") / p.get<int>("N")},
      N{p.get<size_t>("N")},
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

void FD::solve(SimState& state) {
    using namespace blaze;

    // Calculate source term as matrix
    UniformVector<double> ones(N, 1.0);
    auto psi2 = real(state.psis % state.psis);
    DynamicMatrix<double, columnMajor> source =
        dx * dx * expand(psi2 * state.lambda, 1);

    gctrs(dl, d, du, du2, ipiv, source);

    // Solution is only determined up to an additive constant. We fix it by
    // requiring a vanishing mean. This is in accordance with the behavior
    // of Poisson::FFT, where the DC mode is set to 0 explicitly
    double mean = 1.0 / N * sum(source);
    state.V = column(map(source, [&mean](double V) { return V - mean; }), 0);
}

FD::~FD() = default;
}  // namespace Poisson
