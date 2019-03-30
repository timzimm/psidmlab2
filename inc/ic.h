#ifndef __IC__
#define __IC__

#include <fstream>
#include <memory>
#include <string>
#include "common.h"
#include "potential.h"

// Forward Declarations
enum class ICType { Density, Powerspectrum };

class ICGenerator {
   private:
    ICType type;
    int N;
    double dx;
    double L;
    int M;
    double rel_threshold;
    std::string source_name;
    mutable std::ifstream ic_file;
    std::unique_ptr<Potential::Algorithm> potential;
    int data_N;

    void psi_from_rho(SimState& state) const;
    void psi_from_power(SimState& state) const;

   public:
    ICGenerator(const Parameters& param);
    void generate(SimState& state) const;
};
#endif
