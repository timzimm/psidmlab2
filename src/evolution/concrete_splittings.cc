#include "splitting.h"

// Flow Operators
#include "evolution/interaction_external_potential.h"
#include "evolution/kinetic.h"

template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::InteractionExternalPotential>;
template class SRKN<Schroedinger::InteractionExternalPotential,
                    Schroedinger::Kinetic>;

