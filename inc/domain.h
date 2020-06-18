#ifndef __DOMAIN__
#define __DOMAIN__
#include <boost/units/base_units/astronomical/parsec.hpp>
#include "parameters_fwd.h"

// The meaning of each member is best understood by giving a serious of
// examples:
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +             (I) PERIODIC BOUNDARY CONDITIONS                    +
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// (i) Even N=6
// ================================
// x-space
// -------
//  -L/2                       dx = L/N                       L/2-dx
//   =                                                           =
//  xmin     xmin+dx     xmin+2dx    xmin+3dx    xmin+4dx      xmax
//   |          |           |           |           |           |
// [ x0         x1          x2          x3          x4          x5 ]
//   |                                                          |
//   ----- Fields are expected to live on these nodes only ------
//
// k-space
// -------
// -dk*N/2                    dk = 2*pi/L                  dk*N/2 - dk
//   =                                                          =
//  kmin     kmin+dk     kmin+2dk    kmin+3dk    kmin+4dk      kmax
//   |          |           |           |           |           |
// [ k0         k1          k2          k3          k4          k5 ]
//
// (i) Odd N=5
// ================================
// x-space
// -------
//  -L/2                  dx = L/N                L/2-dx
//   =                                              =
//  xmin     xmin+dx     xmin+2dx    xmin+3dx      xmax
//   |          |           |           |           |
// [ x0         x1          x2          x3          x4 ]
//
// k-space
// -------
// NOTE: Truncating integer division (round towards zero)
//
// -dk*N/2               dk = 2*pi/L              dk*N/2
//   =                                              =
//  kmin     kmin+dk     kmin+2dk    kmin+3dk      kmax
//   |          |           |           |           |
// [ k0         k1          k2          k3          k4 ]
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +            (II) HOMOGENEOUS DIRCHLET CONDITIONS                 +
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// (i) Even/Odd N
// ================================
// x-space
// -------
//   dx                       dx = L/(N+1)                     L-dx
//   =                                                          =
//  xmin     xmin+dx     xmin+2dx    xmin+3dx   xmin+(N-2)dx   xmax
//   |          |           |           |           |           |
// [ x0         x1          x2          x3   ...   xN-2       xN-1 ]
//
// k-space
// -------
//
//   dk                         dk = pi/L                      N*dk
//   =                                                          =
//  kmin     kmin+dk     kmin+2dk    kmin+3dk    kmin+(N-2)dk   kmax
//   |          |           |           |           |           |
// [ k0         k1          k2          k3         kN-2         kN-1 ]
//

class Domain {
    // Conversion functions
    double chi_of_x(const boost::units::quantity<
                    boost::units::astronomical::parsec_base_unit::unit_type>
                        x) const;
    double x_of_chi(double chi) const;

   public:
    enum class BoundaryCondition { Periodic, HomogeneousDirichlet };
    const BoundaryCondition bc;
    const bool physical;
    double omega_m0;
    int N;          // No. of NONREDUNDANT spatial points
    double L_phys;  // domain size (Mpc)
    double L;       // domain size (dimless)
    double dx;      // Resolution in x space (dimless)
    double dk;      // Resolution in k space (dimless)
    double m22;     // Mass in units of 1e-22  eV
    double hubble;  // hubble parameter

    // Grid Limits
    double xmin;  // (dimless, nonredundant)
    double xmax;  // (dimless, nonredundant)
    double kmin;  // (dimless, nonredundant)
    double kmax;  // (dimless, nonredundant)
    Domain(const Parameters& p);
};
#endif
