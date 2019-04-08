#ifndef __POISSON_FFT__
#define __POISSON_FFT__

#include <complex>
#include "blaze/math/CompressedMatrix.h"
#include "blaze/math/DiagonalMatrix.h"
#include "blaze/math/DynamicVector.h"
#include "fftw3.h"
#include "interfaces.h"

// Solves Poisson Equation by taking a discrete Fourier transformation of the
// source, transforming its coefficient and taking its inverse FT.
// Assumes periodic boundary conditions.

namespace Poisson {

class FFT : public PotentialMethod::Registrar<FFT> {
    using RCV = blaze::DynamicVector<double>;
    using RRV = blaze::DynamicVector<double, blaze::rowVector>;
    using CCV = blaze::DynamicVector<std::complex<double>>;
    using RDM = blaze::DiagonalMatrix<blaze::CompressedMatrix<double>>;
    size_t N;
    double L;

    CCV fft;
    RDM inv_k_sq;

    fftw_plan forwards;
    fftw_plan backwards;

   public:
    constexpr static double epsilon = 1e-8;

    FFT(const Parameters& p);
    ~FFT();
    void operator()(SimState& state) override;
};

}  // namespace Poisson
#endif
