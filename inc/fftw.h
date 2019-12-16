#ifndef __FFTW__
#define __FFTW__

#include <fftw3.h>
#include <memory>

struct fftw_plan_deleter {
    void operator()(fftw_plan plan) { fftw_destroy_plan(plan); };
};

using fftw_plan_ptr =
    std::unique_ptr<std::remove_pointer_t<fftw_plan>, fftw_plan_deleter>;

#endif
