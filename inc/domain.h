#ifndef __DOMAIN__
#define __DOMAIN__
#include "parameters_fwd.h"

struct Domain {
    const int N;      // No. of spatial points / Fourier basis functions
    const double L;   // domain size (including redundant point at boundary)
    const double dx;  // Resolution in x space
    const double dk;  // Resolution in k space

    // Grid Limits

    // in x-direction we have:
    // [xmin,xmax] = [-L/2, L/2 - dx] (x=L/2 redundant because of periodicity)
    // Note: xmax - xmin != L
    const double xmin;
    const double xmax;
    // in k-direction we have (truncating division):
    // [kmin,kmax] = [-dk*(N/2), dk*(N/2)] if N is odd (idx symmetric around 0)
    // [kmin,kmax] = [-dk*(N/2) + dk, dk*(N/2)] if N is even
    const double kmin;
    const double kmax;
    Domain(const Parameters& p);
};
#endif
