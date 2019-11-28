#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "convolution_functions.h"
#include "fftw3.h"
#include "interfaces.h"

namespace Observable {

class DensityContrast : public ObservableFunctor {
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
    DensityContrast(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(DensityContrast)
};

class PhaseSpaceDistribution : public ObservableFunctor {
    // Abbreviations
    using CCM = blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>;
    using RCM = blaze::DynamicMatrix<double, blaze::columnMajor>;
    using RRM = blaze::DynamicMatrix<double>;
    using RRV = blaze::DynamicVector<double>;

    double sigma_x;  // spatial smoothing scale
    bool husimi;     // compute husimi?
    bool linear;     // do linear convolution?
    int N;           // number of spatial gridpoints
    double dx;       // spatial resolution
    int N_kernel;    // symmetric 5-sigma_x interval in points
    double t_prev;   // timestamp of last cached observable
    convolution_ws<std::complex<double>> ws;
    RCM wigner_f;  // cached wigner
    RRM husimi_f;  // cached husimi
    RRV idx;
    CCM iaf;        // instantaneous autocorrelation function
    fftw_plan c2c;  // IAF -> wigner transform (complex)

    void wigner_distribution(const SimState& state);

    void husimi_distribution(const SimState& state);

   public:
    PhaseSpaceDistribution(const Parameters& p, const Cosmology&);
    ~PhaseSpaceDistribution();
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(PhaseSpaceDistribution)
};

class Potential : public ObservableFunctor {
    std::unique_ptr<Interaction> pot;
    blaze::DynamicVector<double, blaze::columnVector> potential;

   public:
    Potential(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(Potential)
};

class WaveFunction : public ObservableFunctor {
   public:
    WaveFunction(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(WaveFunction)
};

class Energy : public ObservableFunctor {
    const Cosmology& cosmo;
    int N;
    double dx;
    blaze::CompressedMatrix<double> grad;
    blaze::DynamicVector<double> energies;
    blaze::DynamicVector<double, blaze::columnVector> x;

   public:
    Energy(const Parameters& p, const Cosmology& cosmo_);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(Energy)
};

class ParticleFlux : public ObservableFunctor {
    int N;
    double dx;
    blaze::CompressedMatrix<double> grad;
    blaze::DynamicVector<double> flux;

   public:
    ParticleFlux(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(const SimState& state) override;
    REGISTER(ParticleFlux)
};

}  // namespace Observable
#endif
