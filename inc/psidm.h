#ifndef __PSIDM__
#define __PSIDM__
#include <complex>
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

using complex = std::complex<double>;
using c_vector = std::vector<complex>;
using d_vector = std::vector<double>;

struct Parameters {
    std::string filename;  // output filename
    double tau_start;      // initial super conformal time
    double tau_end;        // final super conformal time
    double a_start;        // initial scalefactor
    double a_end;          // final scalefactor
    size_t N;              // Number of spatial points
    size_t M;              // Number of wavefunctions
    double dx;             // spatial resolution
    double L;              // box size
    double mu;             // phase space resolution
};

struct SimState {
    c_vector psis;
    d_vector Vs;
    int n;     // time step number
    double t;  // current time
    double a;  // current scale factor
};

#endif
