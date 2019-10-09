#ifndef __COSMOLOGY__
#define __COSMOLOGY__
#include <boost/units/base_units/astronomical/parsec.hpp>
#include <unordered_map>
#include <vector>

// Forward Declaration
#include "parameters_fwd.h"

// cosmological model "tags"
// Dynamic <=> a = a(tau); omega_m(a) != const
// Artificial <=> tau = tau(a), a = a(tau) given by external file;
enum class CosmoModel { Dynamic, Artificial };

class Cosmology {
   private:
    CosmoModel model;
    double omega_m0;
    double hubble;
    double mu;
    double a_start;
    double a_end;
    double delta_a;
    int A;

    std::vector<double> a_grid;
    std::unordered_map<double, double> tau_a_map;

    double E(const double a) const;
    double dtau_da(const double a) const;

   public:
    Cosmology(const Parameters& p);
    // Time dependent matter density parameter as a function of scalefactor
    double omega_m(double a) const;

    // Linear growth factor
    double Dplus(double a) const;

    // Super conformal time as function of sclaefactor
    double tau_of_a(const double a) const;

    // Inverse of the above
    double a_of_tau(double tau) const;

    // EqualityComparable to the models defined above
    bool operator==(const CosmoModel& model) const;

    // Conversion functions
    double chi_of_x(const boost::units::quantity<
                    boost::units::astronomical::parsec_base_unit::unit_type>
                        x) const;
    double x_of_chi(double chi) const;

    static double z_of_a(const double a) { return 1.0 / a - 1; };
    static double a_of_z(const double z) { return 1.0 / (z + 1); };
};
Parameters& operator<<(Parameters& p, const Cosmology& cosmo);

#endif
