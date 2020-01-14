#pragma once
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

// C++ wrapper to allow capturing lambdas in GSL integration routines
template <typename F>
class gsl_function_pp : public gsl_function {
    const F func;
    static double invoke(double x, void* params) {
        return static_cast<gsl_function_pp*>(params)->func(x);
    }

   public:
    gsl_function_pp(const F& f) : func(f) {
        function = &gsl_function_pp::invoke;
        params = this;
    }
    operator gsl_function*() { return this; }
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
    return gsl_function_pp<F>(func);
}
