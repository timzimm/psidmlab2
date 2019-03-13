#include "common.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iomanip>

Parameters::Parameters(const std::string& filename) {
    namespace pt = boost::property_tree;
    pt::read_ini(filename, tree);

    // Location information
    out_file = tree.get<std::string>("General.output_file");
    ps_file = tree.get<std::string>("General.ps_file");

    // Categorical information
    ic = static_cast<IC>(tree.get<int>("Initial_Conditions.type"));
    cosmo = static_cast<Cosmology>(tree.get<int>("Simulation.cosmology"));
    integrator = static_cast<Integrator>(tree.get<int>("Simulation.algorithm"));

    // Numerical information
    mu = tree.get<double>("Simulation.mu");
    L = tree.get<double>("Simulation.L");
    N = tree.get<size_t>("Simulation.N");
    a_start = tree.get<double>("Simulation.a_start");
    a_end = tree.get<double>("Simulation.a_end");
    dtau = tree.get<double>("Simulation.dtau");
    M = tree.get<size_t>("Initial_Conditions.M");
    ev_thr = tree.get<double>("Initial_Conditions.ev_threshold");

    // inferred information
    dx = L / (N - 1);
    // TODO define a tau conversion
    /* tau_start = tau_of_a(a_start); */
    /* tau_end = tau_of_a(a_end); */
}

std::ostream& operator<<(std::ostream& stream, const Parameters& param) {
    stream << INFOTAG("Simulation Parameters") << std::endl;
    for (const auto& section : param.tree) {
        for (const auto& parameter : section.second) {
            stream << std::left << std::setw(20) << parameter.first << "|"
                   << std::right << std::setw(20) << parameter.second.data()
                   << std::endl;
        }
    }
    stream << std::left;
    return stream;
}

