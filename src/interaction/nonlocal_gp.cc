#include "interaction/nonlocal_gp.h"
#include "fftw3.h"
#include "interfaces.h"
#include "parameters.h"
#include "state.h"

NonlocalGP::NonlocalGP(const Parameters &p, const SimState &state)
    : poisson(p, state), gamma{p["Simulation"]["interaction"]["gamma"]} {}

void NonlocalGP::solve(SimState &state) {
    // compute gravitational potential, i.e. delegate the work.
    poisson.solve(state);
    // add local interaction
    state.V += rho_from(state);
}

void NonlocalGP::solve(blaze::DynamicVector<double> &V,
                       const blaze::DynamicVector<double> &source) {
    // compute gravitational potential
    poisson.solve(V, source);
    // add local interaction
    V += source;
}

