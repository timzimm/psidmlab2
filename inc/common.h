#ifndef __COMMON__
#define __COMMON__
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <complex>
#include <ostream>
#include <string>
#include <vector>

/*This header defines all required simulation parameters as well as what we
 * consider a "simulation state" that gets passed around the code
 * We use structs as we really don't care bout encapsulation at this point. Both
 * objects are merely used to circumvent very long argument lists of the type
 *       void f(vector psi, vectot V, double t, double a, int N,...)
 *  in favor of
 *       voif f(SimState s, Parameters p)
 * We also define some global type abbreviations
 */

// Some convience macros and typedefs no one wants to write out
#define INFOTAG(message) "\033[1;34m[INFO]\033[0m " << message
#define WARNINGTAG(message) "\033[1;33m[WARNING]\033[0m " << message
#define ERRORTAG(message) "\031[1;33m[ERROR]\033[0m " << message

using complex = std::complex<double>;
using c_vector = std::vector<complex>;
using d_vector = std::vector<double>;

struct Parameters {
    enum class IC { Spherical, GRF };
    enum class Cosmology { Static, EDS, LCDM };
    enum class Integrator { USOFFT, USOLW };
    std::string out_file;   // output filename
    std::string ps_file;    // powerspectrum filename
    IC ic;                  // initial condition type
    Cosmology cosmo;        // cosmological model
    Integrator integrator;  // integration algorithm
    double tau_start;       // initial super conformal time
    double tau_end;         // final super conformal time
    double dtau;            // time increment
    double a_start;         // initial scalefactor
    double a_end;           // final scalefactor
    size_t N;               // Number of spatial points
    size_t M;               // Number of wavefunctions
    double ev_thr;          // eigenvalue threshold
    double dx;              // spatial resolution
    double L;               // box size
    double mu;              // phase space resolution
    boost::property_tree::ptree tree;

    Parameters(const std::string& filename);
};

std::ostream& operator<<(std::ostream& stream, const Parameters& param);

struct SimState {
    c_vector psis;
    d_vector Vs;
    int n;       // time step number
    double tau;  // current time
    double a;    // current scale factor
};

#endif
