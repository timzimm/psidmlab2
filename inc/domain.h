#ifndef __DOMAIN__
#define __DOMAIN__
#include "logging.h"
#include "parameters_fwd.h"

#include <array>
#include <vector>

struct Domain {
    using interval_t = std::array<double, 2ul>;

    std::vector<int> Ns;          // degrees of freedom per dimension
    std::vector<interval_t> box;  // box limits as in config (code or physical)
    std::vector<double> box_lengths;  // box lengths (code units)
    std::vector<double> dxs;          // dx's (code units)
    size_t dims;                      // # of dimensions
    size_t N_total;                   // total # of Ns (N_1 * ... * N_dims)

    Domain(const Parameters& p);
};

#endif
