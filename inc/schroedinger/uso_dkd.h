#ifndef __SCHROEDINGER_USO_DKD__
#define __SCHROEDINGER_USO_DKD__

#include <fftw3.h>
#include <complex>
#include <memory>
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "cosmology.h"
#include "interfaces.h"

// Unitary split operator. Propagates psi according to
// psi_(n+1) = exp(-i/2 * V) exp(-i * K) exp(-i/2 * V) psi_(n)
// for separable Hamiltionian H = K + V = -1/2*del^2_x + a(t)*V(psi)

namespace Schroedinger {

class USO_DKD : public SchroedingerMethod::Registrar<USO_DKD> {
    // Abbreviations
    using cmplx = std::complex<double>;

    // Complex row vector
    using CRV = blaze::DynamicVector<cmplx, blaze::rowVector>;
    // Complex column vector
    using CCV = blaze::DynamicVector<cmplx>;
    // Real column vector
    using RCV = blaze::DynamicVector<double>;
    // Complex dense matrix in row-major order
    using CRM = blaze::DynamicMatrix<cmplx>;
    // Complex sparse diagonal matrix in row-major order
    using CDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<cmplx>>;

    // Data members
    Cosmology cosmo;                       // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // Potential method
    size_t N;                              // No. of spatial grid points
    double L;                              // Domain size
    CDM K;  // Kick operator as sparse diagonal matrix
    CDM D;  // Drift operator as sparse diagonal matrix
    RCV wavenumbers;

    fftw_plan forwards;   // in-place forward FFT
    fftw_plan backwards;  // in-place backward FFT

    // Transforms each row of matrix_in according to the passed plan and stores
    // the result in matrix_out
    void transform_matrix(const fftw_plan& plan, CRM& matrix_in,
                          CRM& matrix_out);

    // Kick Operator - returns a blaze expression
    auto kick(const CRM& psis_in_k, const double dt, const double weight);

    // Drift Operator - returns a blaze expression
    auto drift(const CRM& psis_in_x, const RCV& V, const double dt,
               const double t, const double weight);

   public:
    USO_DKD(const Parameters& p);
    ~USO_DKD();
    void operator()(SimState& state) override;
};
}  // namespace Schroedinger
#endif
