#include "cosmology.h"
#include "logging.h"
#include "parameters.h"

#include <boost/math/tools/roots.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/prefixes.hpp>
#include <iostream>
#include <tuple>

using namespace boost::units;
using namespace boost::units::si;
using namespace boost::units::si::constants::codata;
BOOST_UNITS_STATIC_CONSTANT(parsec, astronomical::parsec_base_unit::unit_type);

Cosmology::Cosmology(const Parameters& p)
    : model{static_cast<CosmoModel>(p["model"].get<int>())},
      omega_m0{p["omega_m0"].get<double>()},
      hubble{p["h"].get<double>()},
      mu{p["mu"].get<double>()},
      a_start{a_of_z(p["z_start"].get<double>())},
      a_end{a_start},
      delta_a{0},
      A{p["a_grid_N"].get<int>()},
      a_grid{},
      tau_a_map{} {
    if (model != CosmoModel::Static) {
        // Here the save_at array in the config file holds redshift values. We
        // transform them to scalefactor values and take the largest of them
        // to determine a_end
        auto save_at = p["save_at"].get<std::vector<double>>();
        std::transform(save_at.begin(), save_at.end(), save_at.begin(), a_of_z);
        a_end = *std::max_element(save_at.begin(), save_at.end());
        delta_a = (a_end - a_start) / (A - 1);

        // Super conformal time is defined via its differential
        // Thus, after integrating we still have to fix the integration
        // constant. We do this by defining tau = 0 for the simulation start
        // time. tau increases strictly monotone so our time parameter can
        // be interpreted as a proper time.
        a_grid.resize(A);
        std::cout << INFOTAG("Initialize time lookup table") << std::endl;
        double tau = 0;
        double a = a_start;
        for (int j = 0; j < A; j++) {
            a_grid[j] = a;
            tau_a_map[a] = tau;
            a += delta_a;
            tau += dtau_da(a - 0.5 * delta_a) * delta_a;
        }
    }
}

double Cosmology::omega_m(double a) const {
    if (model == CosmoModel::Static) return omega_m0;
    return omega_m0 / (omega_m0 + (1 - omega_m0) * a * a * a);
}
double Cosmology::Dplus(double a) const {
    double om = omega_m(a);
    return 5.0 * a / 2 * om /
           (pow(om, 4.0 / 7) - (1 - om) + (1 + 0.5 * om) * (1 + (1 - om) / 70));
}

double Cosmology::E(const double a) const {
    return sqrt(omega_m0 / (a * a * a) + (1 - omega_m0));
}

double Cosmology::dtau_da(const double a) const {
    return 1 / (a * a * a * E(a)) * sqrt(1.5 * omega_m0);
}

double Cosmology::tau_of_a(const double a) const {
    // Compute index of closest expansion factor that is strored in the
    // hashtable and smaller than a
    int j;

    // Edge cases extrapolation
    if (a <= a_start)
        return tau_a_map.at(a_grid[0]);
    else if (a >= a_end)
        return tau_a_map.at(a_grid[A - 1]);
    // Linear interpolation
    else
        j = static_cast<int>((a - a_grid[0]) / (a_grid[1] - a_grid[0]));

    // Compute differences to bounding interval
    double lower_a = a_grid[j];
    double upper_a = a_grid[j + 1];
    double lower_delta = a - lower_a;
    double upper_delta = upper_a - a;

    // Compute tau(a) by linear interpolation
    return (lower_delta * tau_a_map.at(upper_a) +
            upper_delta * tau_a_map.at(lower_a)) /
           (lower_delta + upper_delta);
}

double Cosmology::a_of_tau(double tau) const {
    using namespace boost::math::tools;

    if (model == CosmoModel::Static) return a_start;

    // Simulation start
    if (tau <= tau_a_map.at(a_grid[0])) return a_start;

    // Simulation end
    if (tau >= tau_a_map.at(a_grid[A - 1])) return a_end;

    // in between start and end
    auto tau_offset = [&](double a) { return tau_of_a(a) - tau; };
    double a_guess = 0.5 * (a_grid[1] - a_grid[0]);
    double factor = 2;
    bool is_rising = true;
    boost::uintmax_t it = 100;
    int digits = std::numeric_limits<double>::digits - 3;
    eps_tolerance<double> tol(digits);
    auto root =
        bracket_and_solve_root(tau_offset, a_guess, factor, is_rising, tol, it);
    return root.second;
}

bool Cosmology::operator==(const CosmoModel& m) const { return model == m; }

double Cosmology::chi_of_x(
    const quantity<astronomical::parsec_base_unit::unit_type> x) const {
    using phasevolume = derived_dimension<length_base_dimension, 2,
                                          time_base_dimension, -1>::type;

    const quantity<energy> eV = e_over_m_e * m_e * volt;
    const quantity<unit<phasevolume, si::system>> mu_si =
        mu * hbar / (1e-22 * eV) * pow<2>(c);
    const auto H0 = static_cast<quantity<frequency>>(hubble * kilo * meter /
                                                     (second * mega * parsec));

    return 1.0 / root<2>(mu_si) * root<4>(1.5 * pow<2>(H0) * omega_m0) *
           static_cast<quantity<length>>(x);
}

double Cosmology::x_of_chi(double chi) const {
    using phasevolume = derived_dimension<length_base_dimension, 2,
                                          time_base_dimension, -1>::type;

    const quantity<energy> eV = e_over_m_e * m_e * volt;
    const quantity<unit<phasevolume, si::system>> mu_si =
        mu * hbar / (1e-22 * eV) * pow<2>(c);
    const auto H0 = static_cast<quantity<frequency>>(hubble * kilo * meter /
                                                     (second * mega * parsec));

    return (root<2>(mu_si) / root<4>(1.5 * pow<2>(H0) * omega_m0) * chi) /
           static_cast<quantity<length>>(1.0 * parsec);
}

Parameters& operator<<(Parameters& p, const Cosmology& cosmo) {
    const bool convert = p["General"]["physical_units"].get<bool>();
    if (convert) {
        const auto L = p["Simulation"]["L"].get<double>() * parsec;
        p["Simulation"]["L_phys"] = L.value();
        p["Simulation"]["L"] = cosmo.chi_of_x(L);

        for (auto& observable : p["Observables"]) {
            if (auto result = observable.find("sigma_x");
                result != observable.end()) {
                const auto sigma_x = result->get<double>() * parsec;
                observable["sigma_phys"] = sigma_x.value();
                observable["sigma_x"] = cosmo.chi_of_x(sigma_x);
            }
        }
    }
    return p;
}
