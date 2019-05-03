#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "convolution_functions.h"
#include "fftw3.h"
#include "interfaces.h"

namespace Observable {

class DensityContrast : public ObservableFunctor::Registrar<DensityContrast> {
    double sigma_x;  // spatial smoothing scale
    bool husimi;
    bool linear;
    int N;
    double dx;
    int N_kernel;
    double t_prev;
    convolution_ws<double> ws;
    blaze::DynamicVector<double> gaussian_kernel;
    blaze::DynamicVector<double> delta;

   public:
    DensityContrast(const Parameters& p);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
};

class PhaseSpaceDistribution
    : public ObservableFunctor::Registrar<PhaseSpaceDistribution> {
    double sigma_x;  // spatial smoothing scale
    bool husimi;
    bool linear;
    int N;
    double L;
    double dx;
    int N_kernel;
    double t_prev;
    convolution_ws<std::complex<double>> ws;
    blaze::DynamicVector<double> gaussian_kernel;
    blaze::DynamicMatrix<double, blaze::columnMajor> f;

    void wigner_distribution(const SimState& state);

    void husimi_distribution(const SimState& state);

   public:
    PhaseSpaceDistribution(const Parameters& p);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
};

class Potential : public ObservableFunctor::Registrar<Potential> {
   public:
    Potential(const Parameters& p);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
};

class WaveFunction : public ObservableFunctor::Registrar<WaveFunction> {
   public:
    WaveFunction(const Parameters& p);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
};

}  // namespace Observable
#endif
