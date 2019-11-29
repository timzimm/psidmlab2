#include "splitting.h"

// Flow Operators
#include "evolution/interaction_potential_trapezodial.h"
#include "evolution/kinetic.h"

template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::InteractionPotentialTrapezodial>;
template class SRKN<Schroedinger::InteractionPotentialTrapezodial,
                    Schroedinger::Kinetic>;

