#ifndef __IC__
#define __IC__
#include "parameters.h"

#include <cassert>
#include <fstream>
#include <memory>
#include <string>

// Forward Declarations
struct SimState;
class Interaction;
class Cosmology;

enum class ICType {
    ExternalRealImag,
    ExternalModulusPhase,
    Powerspectrum,
    PreviousSimulation
};

class ICGenerator {
   private:
    ICType type;
    int N;
    int data_N;
    double L;
    double dx;
    int seed;
    bool compute_velocity;
    Parameters param;  // This is awfull!
    std::string filename;
    mutable std::ifstream ic_file;

    // init wavefunction ...
    // according to a matter power spectrum provided by file.
    void delta_from_power(SimState& state, const Cosmology& cosmo) const;

    // init wavefunction from file
    void real_imag_from_file(SimState& state) const;
    void modulus_phase_from_file(SimState& state) const;

    // restore wavefunction from previous simulation
    void psi_from_state(SimState& state) const;

   public:
    ICGenerator(const Parameters& param);
    void generate(SimState& state, const Cosmology&) const;

    ~ICGenerator();
};
#endif
