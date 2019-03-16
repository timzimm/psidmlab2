#include "common.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iomanip>
#include "cosmology.h"

Parameters::Parameters(const std::string& filename) {
    namespace pt = boost::property_tree;
    pt::read_ini(filename, tree);

    // Location information
    out_file = tree.get<std::string>("General.output_file");
    ps_file = tree.get<std::string>("General.ps_file");

    // Categorical information
    ic = static_cast<IC>(tree.get<int>("Initial_Conditions.ic_type"));
    cosmo = static_cast<Model>(tree.get<int>("Simulation.cosmology"));
    integrator = static_cast<Integrator>(tree.get<int>("Simulation.algorithm"));

    // Numerical information
    mu = tree.get<double>("Simulation.mu");
    L = tree.get<double>("Simulation.L");
    N = tree.get<size_t>("Simulation.N");
    dtau = tree.get<double>("Simulation.dtau");
    M = tree.get<size_t>("Initial_Conditions.M");
    ev_thr = tree.get<double>("Initial_Conditions.ev_threshold");
    A = tree.get<size_t>("Simulation.A");
    tau_start = 0;
    tau_end = tree.get<double>("Simulation.t_end");
    a_start = Cosmology::a_of_z(tree.get<double>("Simulation.z_start"));

    // inferred information
    dx = L / (N - 1);
    omega_m0 =
        (cosmo == Model::EDS) ? 1 : tree.get<double>("Simulation.omega_m0");
    a_end = (cosmo == Model::Static)
                ? a_start
                : Cosmology::a_of_z(tree.get<double>("Simulation.z_end"));
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

