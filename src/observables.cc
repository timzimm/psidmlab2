#include "observables.h"
#include "blaze/math/UniformMatrix.h"
#include "parameters.h"
#include "state.h"

namespace Observable {

DensityContrast::DensityContrast(const Parameters& p){};

blaze::DynamicMatrix<double, blaze::columnMajor> DensityContrast::compute(
    const SimState& state) {
    // drho(x,t) can directly be constructed from the state.
    auto psi2 = blaze::real(state.psis % state.psis);

    // density constrast delta = |psi|^2 + 1 pointwise
    blaze::UniformVector<double> one(psi2.rows());
    auto delta = psi2 * state.lambda + one;
    return blaze::expand(delta, 1);
}

}  // namespace Observable
