#ifndef __SCHROEDINGER_PCCN__
#define __SCHROEDINGER_PCCN__

#include <complex>
#include <memory>
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "blaze/math/SymmetricMatrix.h"
#include "cosmology.h"
#include "interfaces.h"

// Steps Schoredingers equation forward by applying a Crank Nicolson scheme
// of the form
//
//(psi^{n+1} - psi^n) / dt = 1/2 * (H[psi^n]*psi^n + H[psi^{n+1}]*psi^{n+1})
//
// To circumvent the non-linear system we apply a predictor-corrector
// appproximation of the form:
//
// (psi^{*} - psi^n) / dt = H[psi^n]*psi^n
// (psi^{*} - psi^n) / dt = 1/2 * (H[psi^{n}]*psi^n + H[psi^{*}]*psi^{*})
//
// Using this approximation we don't even have to solve a linear system anymore

namespace Schroedinger {

class PCCN : public SchroedingerMethod::Registrar<PCCN> {
    // complex column vector
    using CCV = blaze::DynamicVector<std::complex<double>>;
    // dense complex matrix
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;
    // sparse symmetric matrix
    using RSM = blaze::SymmetricMatrix<blaze::CompressedMatrix<double>>;

    const Cosmology& cosmo;
    double dx;
    size_t N;
    std::unique_ptr<PotentialMethod> potential;
    RSM K;

   public:
    PCCN(const Parameters& p, const SimState& state, const Cosmology& cosmo_);
    void step(SimState& state) override;
};

}  // namespace Schroedinger
#endif
