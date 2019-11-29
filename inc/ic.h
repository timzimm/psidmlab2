#ifndef __IC__
#define __IC__

#include <cassert>
#include <fstream>
#include <memory>
#include <string>

// Forward Declarations
#include "parameters_fwd.h"
struct SimState;
class Interaction;
class Cosmology;

enum class ICType { ExternalDelta, ExternalPsi, Powerspectrum };

class ICGenerator {
   private:
    ICType type;
    int N;
    int data_N;
    double L;
    double dx;
    int seed;
    double rel_threshold;
    bool compute_velocity;
    mutable std::ifstream ic_file;
    std::unique_ptr<Interaction> pot;

    // init wavefunction ...
    //
    // by loading a file containing delta(x_i)
    void delta_from_file(SimState& state) const;

    // according to a matter power spectrum provided by file.
    void delta_from_power(SimState& state, const Cosmology& cosmo) const;

    // init wavefunction from file
    void psi_from_file(SimState& state) const;

   public:
    ICGenerator(const Parameters& param);
    void generate(SimState& state, const Cosmology&) const;

    // Required to deal with incomplete Interaction type
    ~ICGenerator();
};
#endif
