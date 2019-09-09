#include "splitting.h"

// Flow Operators
#include "schroedinger/cap_poisson_potential.h"
#include "schroedinger/cap_poisson_potential_delta.h"
#include "schroedinger/kinetic.h"
#include "schroedinger/poisson_potential.h"

template class SRKN<Schroedinger::Kinetic, Schroedinger::PoissonPotential>;
template class SRKN<Schroedinger::Kinetic, Schroedinger::CAPPoissonPotential>;
template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::CAPPoissonPotentialDelta>;
