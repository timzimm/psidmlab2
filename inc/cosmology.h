#ifndef __COSMOLOGY__
#define __COSMOLOGY__
#include <unordered_map>
#include <vector>

// Forward Declaration
#include "parameters_fwd.h"

// cosmological model "tags"
// Static <=> a = a_start = const; omega_m(a) = omega_m0
// Dynamic <=> a = a(tau); omega_m(a) != const
enum class CosmoModel { Static, Dynamic };

class Cosmology {
   private:
    CosmoModel model;
    double omega_m0;
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

    // Conversion functions from scalefactor to redshift and vice versa
    static double z_of_a(const double a) { return 1.0 / a - 1; };
    static double a_of_z(const double z) { return 1.0 / (z + 1); };
};

#endif
