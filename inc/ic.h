#ifndef __IC__
#define __IC__

#include <fstream>
#include <string>
#include "common.h"

class ICGenerator {
   private:
    Parameters::IC type;
    int N;
    double dx;
    double L;
    int M;
    double rel_threshold;
    std::string source_name;
    mutable std::ifstream ic_file;
    int data_N;

    void psi_from_rho(SimState& state) const;
    void psi_from_power(SimState& state) const;

   public:
    ICGenerator(const Parameters& param);
    void generate(SimState& state) const;
};
#endif
