#include "domain.h"
#include "parameters.h"

Domain::Domain(const Parameters& p)
    : N{p["Domain"]["N"]},
      L{p["Domain"]["L"]},
      dx(L / N),
      dk(2 * M_PI / L),
      xmin(-L / 2),
      xmax(-xmin - dx),
      // Recall that k is also periodic so we can equally adjust the
      // right boundary:
      //
      // [-dk*(N/2), dk*(N/2)] if N is odd
      // [-dk*(N/2), dk*(N/2) - dk] if N is even
      //
      // We do this to have a "common lower boundary" for both the even and odd
      // case.
      kmin(-dk * (N / 2)),
      kmax((N % 2) ? -kmin : -kmin - dk) {}
