#include <fftw3.h>
#include <cassert>
#include <complex>

class FFT : PoissonSolver<FFT> {
    using namespace blaze;

    DynamicVector<double, rowVector> potential;
    DynamicVector<std::complex<double>> potential_fft;

    DynamicVector<std::complex<double>> source_fft;
    DiagonalMatrix<DynamicMatrix<double>> inv_k_sq;

    double L;
    fftw_plan backwards;

   public:
    static double epsilon = 1e-8;

    FFT(const size_t N, const double L_)
        : potential(N),
          potential_fft(N / 2 - 1),
          source_fft(N / 2 - 1),
          inv_k_sq(N / 2 - 1),
          L(L_) {
        auto diag = diagonal(inv_k_sq);
        // makes the DC constraint manifest
        std::iota(diag.begin(), diag.end(), 0);
        // TODO: Improve!
        diag = -L * L / (4 * M_PI * M_PI * N) * pow(diag, -2);

        backwards = fftw_plan_dft_c2r_1d(
            N, reinterpret_cast<fftw_complex*>(potential_fft.data()), potential,
            FFTW_ESTIMATE);
    }

    ~FFT(){fftw_destroy_plan(backwards)};

    DynamicVector<double, rowVector> solve(
        const DynamicVector<double, rowVector>& source) {
        // Unfortunately, we have to regenerate the transformation plan on each
        // call since the output array location is not fixed. This is still
        // better than than pass by value though.

        // TODO check out fftw_execute_dft_r2c and if alignment requirements are
        // satisfied
        auto forwards = fftw_plan_dft_r2c_1d(
            N, source.data(),
            reinterpret_cast<fftw_complex*>(source_fft.data()), FFTW_ESTIMATE);
        fftw_execute(forwards);
        fftw_destroy_plan(forwards);

        // Check if the DC constrained is satisfied
        assert(std::abs(source_fft[0]) < epsilon);

        // Compute fourier coefficients of the potential
        potential_fft = inv_k_sq * source_fft;

        // Get potential by applying the IFFT
        fftw_execute(backwards);

        // TODO: check if RTO applies. If not, std::move on return and resize on
        // entry
        return potential;
    }
}
