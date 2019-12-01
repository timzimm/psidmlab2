#ifndef __IC__
#define __IC__

// Forward Declarations
#include "parameters_fwd.h"
struct SimState;
class Cosmology;

enum class ICType { ExternalDelta, ExternalPsi, Powerspectrum };

void generate_ic(SimState& state, const Cosmology& cosmo, const Parameters& p);

#endif
