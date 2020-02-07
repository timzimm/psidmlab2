#include "domain.h"
#include "logging.h"
#include "parameters.h"

#include <array>
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
      box_lengths{},
      box_lengths_pc{},
      dxs{},
      dxs_pc{},
      dims{static_cast<int>(Ns.size())},
      N_total{0} {
    using interval_t = std::array<double, 2ul>;
    std::vector<interval_t> box;
    for (auto N : Ns) {
        if (N <= 0) {
            std::cout << ERRORTAG("Ill-defined dofs") << std::endl;
            exit(1);
        }
    }
    for (const auto& interval : p["Domain Properties"]["limits"]) {
        box.emplace_back(interval);
        if (box.back()[1] - box.back()[0] <= 0) {
            std::cout << ERRORTAG("Ill-defined limits") << std::endl;
            exit(1);
        }
        box_lengths.emplace_back(box.back()[1] - box.back()[0]);
        box_lengths_pc.emplace_back(box.back()[1] - box.back()[0]);
    }
    if (box.size() != dims) {
        std::cout << ERRORTAG("Ill-defined box") << std::endl;
        exit(1);
    }

    const double h{p["Cosmology"]["h"]};
    const double omega_m0{p["Cosmology"]["omega_m0"]};
    const double mu{p["Cosmology"]["mu"]};

    // Unit conversion is error prone. We rely on boost:units for
    // that matter to have compile time assurance for correct
    // dimensions.
    using phaseV = derived_dimension<length_base_dimension, 2,
                                     time_base_dimension, -1>::type;
    const quantity<energy> eV = e_over_m_e * m_e * volt;
    const quantity<unit<phaseV, si::system>> mu_si =
        mu * hbar / (1e-22 * eV) * pow<2>(c);
    const auto H0 = static_cast<quantity<frequency>>(h * 100 * kilo * meter /
                                                     (second * mega * parsec));
    const auto timescale = 1.0 / root<2>(1.5 * pow<2>(H0) * omega_m0);

    // Unit Conversion
    if (p["Domain Properties"]["physical_units"].get<bool>()) {
        auto chi_of_x = [&](const double x_val) {
            return 1.0 / root<2>(mu_si * timescale) *
                   static_cast<quantity<length>>(x_val * parsec);
        };
        std::transform(box_lengths_pc.begin(), box_lengths_pc.end(),
                       box_lengths.begin(), chi_of_x);
    } else {
        auto x_of_chi = [&](const double chi) {
            return root<2>(mu_si * timescale) /
                   static_cast<quantity<length>>(1.0 * parsec) * chi;
        };
        std::transform(box_lengths.begin(), box_lengths.end(),
                       box_lengths_pc.begin(), x_of_chi);
    }

    // Grid Point Adjustments
    // Only allow even numbered grids in each dimension
    auto next_even = [](int i) {
        return static_cast<int>(std::round(i / 2.0) * 2.0);
    };
    std::transform(Ns.begin(), Ns.end(), Ns.begin(), next_even);

    // Expand to 3D but set number of grid points to 1 in the new dimensions
    Ns.resize(3, 1);

    // Spatial Resolutions
    for (size_t i = 0; i < dims; ++i) {
        dxs[i] = box_lengths[i] / Ns[i];
        dxs_pc[i] = box_lengths_pc[i] / Ns[i];
    }

    N_total = std::accumulate(Ns.begin(), Ns.end(), 1, std::multiplies<>());
}

