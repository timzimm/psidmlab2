#include "fftw.h"

#ifdef PSIDMLAB_SMP
#include <omp.h>
#include <algorithm>
#endif

// Simple wrappers for all apearing FFT types.
// Complex-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft(int N, const fftw_complex* vin,
                                 const fftw_complex* vout, int dir,
                                 unsigned flags) {
    auto in = const_cast<fftw_complex*>(vin);
    auto out = const_cast<fftw_complex*>(vout);

#ifdef PSIDMLAB_SMP
    const int max_N_per_thread = 1 << 15;
    const int nthreads = std::min(N / max_N_per_thread, omp_get_max_threads());
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_1d(N, in, out, dir, flags));
}
// Real-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft_r2c(int N, const double* vin,
                                     const fftw_complex* vout, unsigned flags) {
    auto in = const_cast<double*>(vin);
    auto out = const_cast<fftw_complex*>(vout);

#ifdef PSIDMLAB_SMP
    const int max_N_per_thread = 1 << 15;
    const int nthreads = std::min(N / max_N_per_thread, omp_get_max_threads());
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_r2c_1d(N, in, out, flags));
}

// Complex-to-Real FFT
fftw_plan_ptr make_fftw_plan_dft_c2r(int N, const fftw_complex* vin,
                                     const double* vout, unsigned flags) {
    auto in = const_cast<fftw_complex*>(vin);
    auto out = const_cast<double*>(vout);

#ifdef PSIDMLAB_SMP
    const int max_N_per_thread = 1 << 15;
    const int nthreads = std::min(N / max_N_per_thread, omp_get_max_threads());
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_dft_c2r_1d(N, in, out, flags));
}

// Real-to-Real transform
fftw_plan_ptr make_fftw_plan_r2r_1d(int n, const double* vin,
                                    const double* vout, fftw_r2r_kind kind,
                                    unsigned flags) {
    auto in = const_cast<double*>(vin);
    auto out = const_cast<double*>(vout);

#ifdef PSIDMLAB_SMP
    const int max_N_per_thread = 1 << 15;
    const int nthreads = std::min(n / max_N_per_thread, omp_get_max_threads());
    fftw_plan_with_nthreads(nthreads);
#endif

    return fftw_plan_ptr(fftw_plan_r2r_1d(n, in, out, kind, flags));
}
