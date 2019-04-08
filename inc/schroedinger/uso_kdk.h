#ifndef __SCHROEDINGER_USO_KDK__
#define __SCHROEDINGER_USO_KDK__

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
// psi_(n+1) = exp(-i/2 * K) exp(-i * V) exp(-i/2 * K) psi_(n)
// for separable Hamiltionian H = K + V = -1/2*del^2_x + a(t)*V(psi)
namespace Schroedinger {

class USO_KDK : public SchroedingerMethod::Registrar<USO_KDK> {
    // Abbreviations
    using cmplx = std::complex<double>;

    // Complex row vector
    using CRV = blaze::DynamicVector<cmplx, blaze::rowVector>;
    // Complex column vector
    using CCV = blaze::DynamicVector<cmplx>;
    // Real column vector
    using RCV = blaze::DynamicVector<double>;
    // Complex dense matrix in column-major order
    using CCM = blaze::DynamicMatrix<cmplx, blaze::columnMajor>;
    // Complex sparse diagonal matrix in row-major order
    using CDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<cmplx>>;

    // Data members
    Cosmology cosmo;                       // cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    size_t N;                              // No. of spatial points
    double L;                              // size of domain
    bool firstStep;                        // is firstStep?
    CDM K;  // Kick operator as sparse diagonal matrix
    CDM D;  // Drift operator as sparse diagonal matrix
    RCV wavenumbers;

    fftw_plan forwards;      // in-place forward FFT
    fftw_plan backwards;     // in-place backward FFT
    fftw_plan forwards_op;   // out-place forward FFT
    fftw_plan backwards_op;  // out-place backward FFT

    // Buffer to store k representation of psi for KDK
    CCM psis_cached;

    // Transforms each row of matrix_in according to the passed plan and stores
    // the result in matrix_out
    void transform_matrix(const fftw_plan& plan, CCM& matrix_in,
                          CCM& matrix_out);

    // Kick Operator - returns a blaze expression
    auto kick(const CCM& psis_in_k, const double dt, const double weight);

    // Drift Operator - returns a blaze expression
    auto drift(const CCM& psis_in_x, const RCV& V, const double dt,
               const double t, const double weight);

   public:
    USO_KDK(const Parameters& p);
    ~USO_KDK();
    void operator()(SimState& state) override;
};
}  // namespace Schroedinger
#endif
