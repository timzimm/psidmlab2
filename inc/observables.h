#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "convolution_functions.h"
#include "fftw3.h"
#include "interfaces.h"

namespace Observable {

class DensityContrast : public ObservableFunctor::Registrar<DensityContrast> {
    double sigma_x;  // spatial smoothing scale
    bool husimi;
    int N;
    double dx;
    int N_kernel;
    convolution_ws ws;
    blaze::DynamicVector<double> gaussian_kernel;

   public:
    DensityContrast(const Parameters& p);
    blaze::DynamicMatrix<double, blaze::columnMajor> compute(
        const SimState& state) override;
};

}  // namespace Observable
#endif
