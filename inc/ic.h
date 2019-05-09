#ifndef __IC__
#define __IC__

#include <cassert>
#include <fstream>
#include <memory>
#include <string>

// Forward Declarations
#include "parameters_fwd.h"
struct SimState;
class PotentialMethod;

enum class ICType { External, Powerspectrum, Density };

class ICGenerator {
   private:
    ICType type;
    int N;
    int data_N;
    double L;
    double dx;
    int M;
    double rel_threshold;
    mutable std::ifstream ic_file;
    mutable std::ifstream pot_file;
    std::unique_ptr<PotentialMethod> potential;

    // init wavefunction matrix ...
    //
    // by loading a file containing Im(psi_i) and Re(psi_i) for all i.
    void psi_from_file(SimState& state) const;

    // by solving the associated operator eigenvalue problem
    // for f(x,t) = rho(x) * delta(0) (cold initial conditions)
    // rho is provided by a file.
    void psi_from_rho(SimState& state) const;

    // according to a matter power spectrum provided by file.
    void psi_from_power(SimState& state) const;

   public:
    ICGenerator(const Parameters& param);
    void generate(SimState& state, Parameters& param) const;

    // Required to deal with incomplete PotentialMethod type
    ~ICGenerator();
};
#endif
