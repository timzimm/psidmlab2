#ifndef __COSMOLOGY__
#define __COSMOLOGY__
#include <unordered_map>
#include <vector>

// Forward Declaration
#include "parameters_fwd.h"

enum class CosmoModel { Static, EDS, LCDM };

class Cosmology {
   private:
    double a_start;
    double a_end;
    double A;
    double delta_a;
    std::vector<double> a_grid;
    CosmoModel model;
    double omega_m0;
    std::unordered_map<double, double> tau_a_map;

    double E(const double a) const;
    double dtau_da(const double a) const;

   public:
    Cosmology(const Parameters& p);
    double tau_of_a(const double a) const;
    double a_of_tau(double tau) const;
    static double z_of_a(const double a);
    static double a_of_z(const double z);
};

#endif
