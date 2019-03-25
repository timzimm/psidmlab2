#ifndef __COMMON__
#define __COMMON__
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <boost/property_tree/ptree.hpp>
#include <complex>
#include <ostream>
#include <string>

/*This header defines all required simulation parameters as well as what we
 * consider a "simulation state" that gets passed around the code
 * We use structs as we really don't care bout encapsulation at this point. Both
 * objects are merely used to circumvent very long argument lists of the type
 *       void f(vector psi, vectot V, double t, double a, int N,...)
 *  in favor of
 *       voif f(SimState s, Parameters p)
 */

// Some convience macros and typedefs no one wants to write out
#define INFOTAG(message) "\033[1;34m[INFO]\033[0m " << message
#define WARNINGTAG(message) "\033[1;33m[WARNING]\033[0m " << message
#define ERRORTAG(message) "\031[1;33m[ERROR]\033[0m " << message

// Forward declarations
namespace Poisson {
enum class Type;
}
enum class ICType;
enum class CosmoModel;

struct Parameters {
    // TODO move enums out here!
    enum class Integrator { USOFFT, USOLW };
    std::string out_file;        // output filename
    std::string ic_source_file;  // powerspectrum filename
    ICType ic;                   // initial condition type
    Poisson::Type psolver;       // Poisson solver type
    CosmoModel cosmo;            // cosmological model
    Integrator integrator;       // integration algorithm
    double tau_start;            // initial super conformal time
    double tau_end;              // final super conformal time
    double dtau;                 // time increment
    double a_start;              // initial scalefactor
    double a_end;                // final scalefactor
    size_t N;                    // Number of spatial points
    size_t M;                    // Number of wavefunctions
    size_t A;                    // Number of points of the scalefactor grid
    double ev_thr;               // eigenvalue threshold
    double dx;                   // spatial resolution
    double L;                    // box size
    double mu;                   // phase space resolution
    double omega_m0;             // matter density parameter
    boost::property_tree::ptree tree;

    Parameters(const std::string& filename);
};

std::ostream& operator<<(std::ostream& stream, const Parameters& param);

struct SimState {
    int n;       // time step number
    double tau;  // current time
    double a;    // current scale factor
    blaze::DynamicVector<double, blaze::rowVector> V;

    // state = sum_i lambda_i * |psi_i><psi_i|
    int M;
    blaze::DynamicMatrix<std::complex<double>> psis;
    blaze::DynamicVector<double, blaze::rowVector> rho_bg;
    blaze::DiagonalMatrix<blaze::DynamicMatrix<double>> lambda;

    SimState(const Parameters& param);
};

#endif
