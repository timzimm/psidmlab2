#include "domain.h"
#include "boost/units/systems/detail/constants.hpp"
#include "logging.h"
#include "parameters.h"

#include <boost/math/tools/roots.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/prefixes.hpp>

using namespace boost::units;
using namespace boost::units::si;
using namespace boost::units::si::constants::codata;
BOOST_UNITS_STATIC_CONSTANT(parsec, astronomical::parsec_base_unit::unit_type);

Domain::Domain(const Parameters& p)
    : bc{static_cast<BoundaryCondition>(p["Domain"]["boundary_type"])},
      physical{p["Domain"]["physical_units"]},
      hubble(0),
      m22(0),
      omega_m0(0),
      N{p["Domain"]["N"]},
      L_phys{p["Domain"]["L"]},
      L(L_phys) {
    if (L <= 0) {
        std::cerr << ERRORTAG("Box has negative length") << std::endl;
        exit(1);
    }
    if (N <= 0) {
        std::cerr << ERRORTAG("Box has negative number of grid points")
                  << std::endl;
        exit(1);
    }
    if (physical) {
        if (m22 <= 0) {
            std::cerr << ERRORTAG("mass parameter is negtive") << std::endl;
            exit(1);
        }
        if (hubble <= 0) {
            std::cerr << ERRORTAG("Box has negative number of grid points")
                      << std::endl;
            exit(1);
        }
        hubble = p["Domain"]["h"];
        m22 = p["Domain"]["m22"];
        omega_m0 = p["Cosmology"]["omega_m0"];
        L = chi_of_x(L_phys * 1e6 * parsec);
    }
    if (bc == BoundaryCondition::Periodic) {
        dx = L / N;
        dk = 2 * M_PI / L;
        xmin = -L / 2;
        xmax = -xmin - dx;
        kmin = -dk * (N / 2);
        kmax = (N % 2) ? -kmin : -kmin - dk;
    } else if (bc == BoundaryCondition::HomogeneousDirichlet) {
        dx = L / (N + 1);
        dk = M_PI / L;
        xmin = dx;
        xmax = L - dx;
        kmin = dk;
        kmax = N * dk;
    }
}

double Domain::chi_of_x(
    const quantity<astronomical::parsec_base_unit::unit_type> x) const {
    using phasevolume = derived_dimension<length_base_dimension, 2,
                                          time_base_dimension, -1>::type;

    const quantity<energy> eV = e_over_m_e * m_e * volt;
    const quantity<unit<phasevolume, si::system>> mu_si =
        hbar / (m22 * 1e-22 * eV) * pow<2>(c);
    const auto H0 = static_cast<quantity<frequency>>(
        hubble * 100 * kilo * meter / (second * mega * parsec));
    const auto alpha = 1.0 / root<2>(1.5 * pow<2>(H0) * omega_m0);

    return 1.0 / root<2>(mu_si * alpha) * static_cast<quantity<length>>(x);
}
