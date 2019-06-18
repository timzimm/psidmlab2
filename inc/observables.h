#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "convolution_functions.h"
#include "fftw3.h"
#include "interfaces.h"

namespace Observable {

class DensityContrast : public ObservableFunctor::Registrar<DensityContrast> {
    const double sigma_x;       // spatial smoothing scale
    const bool husimi;          // compute husimi?
    const bool linear;          // do linear convolution?
    const int N;                // number of spatial gridpoints
    const double dx;            // spatial resolution
    const int N_kernel;         // symmetric 5-sigma_x interval in points
    double t_prev;              // timestamp of last cached observable
    convolution_ws<double> ws;  // holds husimi delta
    blaze::DynamicVector<double> delta;  // holds wigner delta

   public:
    DensityContrast(const Parameters& p);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
};

class PhaseSpaceDistribution
    : public ObservableFunctor::Registrar<PhaseSpaceDistribution> {
    // Abbreviations
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;
    using RCM = blaze::DynamicMatrix<double, blaze::columnMajor>;
    using RRM = blaze::DynamicMatrix<double>;

    double sigma_x;  // spatial smoothing scale
    bool husimi;     // compute husimi?
    bool linear;     // do linear convolution?
    int N;           // number of spatial gridpoints
    double dx;       // spatial resolution
    int N_kernel;    // symmetric 5-sigma_x interval in points
    double t_prev;   // timestamp of last cached observable
    convolution_ws<std::complex<double>> ws;
    RCM wigner_f;   // cached wigner
    RRM husimi_f;   // cached husimi
    CCM iaf;        // instantaneous autocorrelation function
    fftw_plan c2c;  // IAF -> wigner transform (complex)

    void wigner_distribution(const SimState& state);

    void husimi_distribution(const SimState& state);

   public:
    PhaseSpaceDistribution(const Parameters& p);
    ~PhaseSpaceDistribution();
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
