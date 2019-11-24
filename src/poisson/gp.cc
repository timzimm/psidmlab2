#include "poisson/gp.h"
#include "parameters.h"
#include "state.h"

#include <cassert>

namespace Poisson {
GP::GP(const Parameters& p) {}
GP::~GP() = default;

void GP::solve(SimState& state) { solve(state.V, delta_from(state)); }

void GP::solve(blaze::DynamicVector<double, blaze::columnVector>& V,
               const blaze::DynamicVector<double, blaze::columnVector>& s) {
    V = std::move(s);
}
}  // namespace Poisson
