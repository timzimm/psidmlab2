#ifndef __COSMOLOGY__
#define __COSMOLOGY__
#include <unordered_map>
#include <vector>

// Forward Declaration
#include "parameters_fwd.h"

enum class CosmoModel { Static, Dynamic };

class Cosmology {
   private:
    const double a_start;
    const double a_end;
    const int A;
    const double delta_a;
    std::vector<double> a_grid;
    const CosmoModel model;
    const double omega_m0;
    std::unordered_map<double, double> tau_a_map;

    double E(const double a) const;
    double dtau_da(const double a) const;

   public:
    Cosmology(const Parameters& p);
    double tau_of_a(const double a) const;
    double a_of_tau(double tau) const;
    bool operator==(const CosmoModel& model) const;
    static double z_of_a(const double a);
    static double a_of_z(const double z);
};

#endif
