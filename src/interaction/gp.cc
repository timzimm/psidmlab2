#include "interaction/gp.h"
#include "parameters.h"
#include "state.h"

GP::GP(const Parameters& p) {}
GP::~GP() = default;

void GP::solve(SimState& state) {
    // No copy
    state.V = delta_from(state);
}

void GP::solve(blaze::DynamicVector<double>& V,
               const blaze::DynamicVector<double>& s) {
    // TODO: Move from const is copy !!!
    V = std::move(s);
}
