#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include "fftw3.h"
#include "interfaces.h"

#include <blaze/math/CompressedMatrix.h>

// Solves Poisson Equation by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.
// Assumes periodic boundary conditions.

namespace Poisson {

class FFT : public PotentialMethod {
    using RCV = blaze::DynamicVector<double>;
    using CCV = blaze::DynamicVector<std::complex<double>>;
    using RDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<double>>;
    const size_t N;
    const double L;
    const double epsilon;

    CCV fft;
    RCV source;
    RDM kernel;

    fftw_plan forwards;
    fftw_plan backwards;

   public:
    FFT(const Parameters& p);
    ~FFT();
    void solve(SimState& state) override;
    void solve(blaze::DynamicVector<double, blaze::columnVector>& V,
               const blaze::DynamicVector<double, blaze::columnVector>& source)
        override;
    REGISTER(FFT)
};

}  // namespace Poisson
#endif
