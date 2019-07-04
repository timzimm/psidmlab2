#include "cosmology.h"
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <tuple>
#include "logging.h"
#include "parameters.h"

Cosmology::Cosmology(const Parameters& p)
    : a_start{a_of_z(p["Simulation"]["z_start"].get<double>())},
      a_end{a_of_z(p["Simulation"]["z_end"].get<double>())},
      A{p["Simulation"]["a_grid_N"].get<int>()},
      delta_a{(a_end - a_start) / (A - 1)},
      a_grid{},
      model{static_cast<CosmoModel>(p["Cosmology"]["model"].get<int>())},
      omega_m0{p["Cosmology"]["omega_m0"].get<double>()},
      tau_a_map{} {
    // Super conformal time is defined via its differential
    // Thus, after integrating we still have to fix the integration constant.
    // We do this by defining tau = 0 for the simulation start time.
    // tau increases strictly monotone so our time parameter can be interpreted
    // as a proper time.

    if (model != CosmoModel::Static) {
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
