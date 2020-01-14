#ifndef __DOMAIN__
#define __DOMAIN__
#include "logging.h"
#include "parameters_fwd.h"

#include <array>
#include <vector>

struct Domain {
    std::vector<int> Ns;  // degrees of freedom per dimension. Always 3D.
    std::vector<double> box_lengths;     // box lengths (code units)
    std::vector<double> box_lengths_pc;  // box lengths (pc)
    std::vector<double> dxs;             // dx's (code units)
    std::vector<double> dxs_pc;          // dx's (pc)
    int dims;                            // # of dimensions
    size_t N_total;                      // total # of Ns (N_1 * ... * N_dims)

    Domain(const Parameters& p);
};

#endif
