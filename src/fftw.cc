#include "fftw.h"
#ifdef PSIMDLAB_SMP
#include <omp.h>
#endif

// Simple wrappers for all apearing FFT types.
// Complex-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft(int N, const fftw_complex* vin,
                                 const fftw_complex* vout, int dir,
                                 unsigned flags) {
    auto in = const_cast<fftw_complex*>(vin);
    auto out = const_cast<fftw_complex*>(vout);

#ifdef PSIMDLAB_SMP
    const int max_N_per_thread = 1 << 14;
    const int nthreads =
        std::min({N / max_N_per_thread, omp_get_max_threads()});
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_1d(N, in, out, dir, flags));
}
// Real-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft_r2c(int N, const double* vin,
                                     const fftw_complex* vout, unsigned flags) {
    auto in = const_cast<double*>(vin);
    auto out = const_cast<fftw_complex*>(vout);

#ifdef PSIMDLAB_SMP
    const int max_N_per_thread = 1 << 14;
    const int nthreads =
        std::min({N / max_N_per_thread, omp_get_max_threads()});
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_r2c_1d(N, in, out, flags));
}

// Complex-to-Real FFT
fftw_plan_ptr make_fftw_plan_dft_c2r(int N, const fftw_complex* vin,
                                     const double* vout, unsigned flags) {
    auto in = const_cast<fftw_complex*>(vin);
    auto out = const_cast<double*>(vout);

#ifdef PSIMDLAB_SMP
    const int max_N_per_thread = 1 << 14;
    const int nthreads =
        std::min({N / max_N_per_thread, omp_get_max_threads()});
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_c2r_1d(N, in, out, flags));
}
