#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "fftw3.h"
#include "interfaces.h"

namespace Observable {

/* class enum QuasiProbability { Wigner, Husimi }; */

class DensityContrast : public ObservableFunctor::Registrar<DensityContrast> {
   public:
    DensityContrast(const Parameters& p);
    blaze::DynamicMatrix<double, blaze::columnMajor> compute(
        const SimState& state) override;
};

}  // namespace Observable
#endif
