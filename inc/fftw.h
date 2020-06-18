#ifndef __FFTW__
#define __FFTW__

#include <fftw3.h>
#include <memory>

struct fftw_plan_deleter {
    void operator()(fftw_plan plan) { fftw_destroy_plan(plan); };
};

using fftw_plan_ptr =
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, fftw_plan_deleter>;

// Simple wrappers for all apearing FFT types.

// Complex-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft(int N, const fftw_complex* vin,
                                 const fftw_complex* vout, int dir,
                                 unsigned flags);
// Real-to-Complex FFT
fftw_plan_ptr make_fftw_plan_dft_r2c(int N, const double* vin,
                                     const fftw_complex* vout, unsigned flags);

// Complex-to-Real FFT
fftw_plan_ptr make_fftw_plan_dft_c2r(int N, const fftw_complex* vin,
                                     const double* vout, unsigned flags);
// Real-to-Real transform
fftw_plan_ptr make_fftw_plan_r2r_1d(int n, const double* vin,
                                    const double* vout, fftw_r2r_kind kind,
                                    unsigned flags);

// Real-to-Real transform (advanced interface)
fftw_plan_ptr make_fftw_plan_many_r2r(int rank, const int* n, int howmany,
                                      const double* in, const int* inembed,
                                      int istride, int idist, const double* out,
                                      const int* onembed, int ostride,
                                      int odist, const fftw_r2r_kind* kind,
                                      unsigned flags);

#endif
