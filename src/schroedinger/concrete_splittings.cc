#include "splitting.h"

// Flow Operators
#include "schroedinger/interaction_external_potential.h"
#include "schroedinger/kinetic.h"

template class SRKN<Schroedinger::Kinetic,
                    Schroedinger::InteractionExternalPotential>;
template class SRKN<Schroedinger::InteractionExternalPotential,
                    Schroedinger::Kinetic>;

