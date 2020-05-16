#include "splitting.h"

// Time Evolution Operators
#include "evolution/interaction_potential.h"
#include "evolution/interaction_potential_cap_trapezodial.h"
#include "evolution/interaction_potential_magnus.h"
#include "evolution/kinetic.h"

template class PRK<Schroedinger::InteractionPotential, Schroedinger::Kinetic>;
template class PRK<Schroedinger::InteractionPotentialMagnus,
                   Schroedinger::Kinetic>;
template class PRK<Schroedinger::InteractionPotentialCAPTrapezodial,
                   Schroedinger::Kinetic>;
