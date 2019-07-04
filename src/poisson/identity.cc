#include "poisson/identity.h"
#include "parameters.h"
#include "state.h"

namespace Poisson {

Identity::Identity(const Parameters& p) {}

void Identity::solve(SimState& state) {}
void Identity::solve(
    blaze::DynamicVector<double, blaze::columnVector>& V,
    const blaze::DynamicVector<double, blaze::columnVector>& source) {}

}  // namespace Poisson
