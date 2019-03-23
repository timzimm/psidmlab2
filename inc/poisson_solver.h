#ifndef __POISSON_SOLVER_INTERFACE__
#define __POISSON_SOLVER_INTERFACE__
#include <blaze/math/DynamicVector.h>
#include "poisson_solvers/fft_solver.h"

// Poisson Solvers are implemented using static polymorphism via the CRTP
// technique. The class template below defines the interface all dervied solvers
// are required to implement.

template <typename DerivedSolver>
class PoissonSolver {
    using namespace blaze;
    DerivedSolver& dispatch() { return *static_cast<DerivedSolver*>(this); }

   public:
    // Functions all PoissonSolvers must have in common
    DynamicVector<double, rowVector> solve(
        const DynamicVector<double, rowVector>& source) {
        return dispatch().solve();
    }
};

#endif
