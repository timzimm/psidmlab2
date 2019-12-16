#include "splitting.h"

// Flow Operators
#include "evolution/cap_poisson_potential_delta.h"
#include "evolution/interaction_potential_trapezodial.h"
#include "evolution/kinetic.h"

template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::InteractionPotentialTrapezodial>;
template class SRKN<Schroedinger::InteractionPotentialTrapezodial,
                    Schroedinger::Kinetic>;
template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::CAPPoissonPotentialDelta>;
template class SRKN<Schroedinger::CAPPoissonPotentialDelta,
                    Schroedinger::Kinetic>;

