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
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // potential method
    int N;                                 // No. of spatial points
    double L;                              // size of domain
    RCV k_squared;

    fftw_plan forwards;      // in-place forward FFT
    fftw_plan backwards;     // in-place backward FFT
    fftw_plan forwards_op;   // out-place forward FFT
    fftw_plan backwards_op;  // out-place backward FFT

    // Buffer to store k representation of psi for KDK
    CCM psis_cached;

   public:
    USO_KDK(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    ~USO_KDK();
    void step(SimState& state) override;
};
}  // namespace Schroedinger
#endif
