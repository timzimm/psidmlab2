#include "domain.h"
#include "boost/units/systems/detail/constants.hpp"
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
    : physical{p["Domain"]["physical_units"]},
      hubble(0),
      m22(0),
      omega_m0(0),
      N{p["Domain"]["N"]},
      L_phys{p["Domain"]["L"]},
      L(L_phys) {
    if (physical) {
        hubble = p["Domain"]["h"];
        m22 = p["Domain"]["m22"];
        omega_m0 = p["Cosmology"]["omega_m0"];
        L = chi_of_x(L_phys * 1e6 * parsec);
    }
    dx = L / N;
    dk = 2 * M_PI / L;
    xmin = -L / 2;
    xmax = -xmin - dx;
    // Recall that k is also periodic so we can equally adjust the
    // right boundary:
    //
    // [-dk*(N/2), dk*(N/2)] if N is odd
    // [-dk*(N/2), dk*(N/2) - dk] if N is even
    //
    // We do this to have a "common lower boundary" for both the even and odd
    // case.
    kmin = -dk * (N / 2);
    kmax = (N % 2) ? -kmin : -kmin - dk;
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
