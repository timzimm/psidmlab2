#ifndef __DOMAIN__
#define __DOMAIN__
#include <boost/units/base_units/astronomical/parsec.hpp>
#include "parameters_fwd.h"

class Domain {
    // Conversion functions
    double chi_of_x(const boost::units::quantity<
                    boost::units::astronomical::parsec_base_unit::unit_type>
                        x) const;
    double x_of_chi(double chi) const;

   public:
    const bool physical;
    double omega_m0;
    int N;          // No. of spatial points / Fourier basis functions
    double L_phys;  // domain size (Mpc)
    double L;       // domain size (dimless)
    double dx;      // Resolution in x space (dimless)
    double dk;      // Resolution in k space (dimless)
    double m22;     // Mass in units of 1e-22  eV
    double hubble;  // hubble parameter

    // Grid Limits

    // in x-direction we have:
    // [xmin,xmax] = [-L/2, L/2 - dx] (x=L/2 redundant because of periodicity)
    // Note: xmax - xmin != L
    double xmin;  // (dimless)
    double xmax;  // (dimless)
    // in k-direction we have (truncating division):
    // [kmin,kmax] = [-dk*(N/2), dk*(N/2)] if N is odd (idx symmetric around 0)
    // [kmin,kmax] = [-dk*(N/2) + dk, dk*(N/2)] if N is even
    double kmin;
    double kmax;
    Domain(const Parameters& p);
};
#endif
