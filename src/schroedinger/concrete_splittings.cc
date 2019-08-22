#include "splitting.h"

// Flow Operators
#include "schroedinger/kinetic.h"
#include "schroedinger/poisson_potential.h"

template class StrangSplitting<Schroedinger::Kinetic,
                               Schroedinger::PoissonPotential>;
/* template class StrangSplitting<PMLKinetic, PoissonPotential>; */
