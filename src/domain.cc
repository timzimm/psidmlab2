#include "domain.h"
#include "parameters.h"

#include <boost/units/base_units/astronomical/parsec.hpp>
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
    : Ns{p["Domain Properties"]["dof"]},
      box{},
      box_lengths{},
      dxs{},
      dims{Ns.size()},
      N_total{0} {
    for (const auto& interval : p["Domain Properties"]["limits"]) {
        box.emplace_back(interval);
        box_lengths.emplace_back(abs(box.back()[1] - box.back()[0]));
    }
    if (box.size() != dims) {
        std::cout << ERRORTAG("Ill-defined box") << std::endl;
        exit(1);
    }
    if (p["Domain Properties"]["physical_units"].get<bool>()) {
        const double h{p["Cosmology"]["h"]};
        const double omega_m0{p["Cosmology"]["omega_m0"]};
        const double mu{p["Cosmology"]["mu"]};

        // Unit conversion is error prone. We rely on boost:units for that
        // matter to have compile time assurance for correct dimensions.
        using phasevolume = derived_dimension<length_base_dimension, 2,
                                              time_base_dimension, -1>::type;

        auto chi_of_x =
            [&](const quantity<astronomical::parsec_base_unit::unit_type> x) {
                const quantity<energy> eV = e_over_m_e * m_e * volt;
                const quantity<unit<phasevolume, si::system>> mu_si =
                    mu * hbar / (1e-22 * eV) * pow<2>(c);
                const auto H0 = static_cast<quantity<frequency>>(
                    h * 100 * kilo * meter / (second * mega * parsec));
                const auto timescale =
                    1.0 / root<2>(1.5 * pow<2>(H0) * omega_m0);

                return 1.0 / root<2>(mu_si * timescale) *
                       static_cast<quantity<length>>(x);
            };
        std::transform(box_lengths.begin(), box_lengths.end(),
                       box_lengths.begin(), chi_of_x);
    }

    for (size_t i = 0; i < dims; ++i) dxs[i] = box_lengths[i] / Ns[i];

    N_total = std::accumulate(Ns.begin(), Ns.end(), 1, std::multiplies<>());
}

