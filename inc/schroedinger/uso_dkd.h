#ifndef __SCHROEDINGER_USO_DKD__
#define __SCHROEDINGER_USO_DKD__

#include <fftw3.h>
#include <memory>

#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
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
    // Complex dense matrix in column-major order
    using CCM = blaze::DynamicMatrix<cmplx, blaze::columnMajor>;
    // Complex sparse diagonal matrix in row-major order
    using CDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<cmplx>>;

    // Data members
    const Cosmology& cosmo;                // Cosmological model for a(t)
    std::unique_ptr<PotentialMethod> pot;  // Potential method
    int N;                                 // No. of spatial grid points
    double L;                              // Domain size
    RCV k_squared;

    fftw_plan forwards;   // in-place forward FFT
    fftw_plan backwards;  // in-place backward FFT

    void step_internal(SimState& state, const double dt);

   public:
    USO_DKD(const Parameters& p, const SimState& state,
            const Cosmology& cosmo_);
    ~USO_DKD();
    void step(SimState& state) override;
};
}  // namespace Schroedinger
#endif
