#ifndef __OBSERVABLES__
#define __OBSERVABLES__

#include "convolution_functions.h"
#include "cosmology.h"
#include "fftw.h"
#include "interfaces.h"
#include "parameters.h"

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
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(DensityContrast)
};

class PhaseSpaceDistribution : public ObservableFunctor {
    using interval_t = std::array<double, 2ul>;
    using patch_t = std::array<interval_t, 2ul>;

    const double sigma_x;  // spatial smoothing scale
    const bool husimi;     // compute husimi?
    const bool linear;     // do linear convolution?
    const patch_t patch;
    int N;               // number of spatial gridpoints
    const double L;      // spatial resolution
    const int N_kernel;  // symmetric 5-sigma_x interval in points
    double t_prev;       // timestamp of last cached observable
    convolution_ws<std::complex<double>> ws;
    blaze::DynamicMatrix<double, blaze::columnMajor> wigner_f;  // cached wigner
    blaze::DynamicMatrix<double> husimi_f;                      // cached husimi
    blaze::DynamicVector<double> idx;
    blaze::DynamicVector<double> idk;
    blaze::DynamicMatrix<std::complex<double>, blaze::columnMajor>
        iaf;            // instantaneous autocorrelation function
    fftw_plan_ptr c2c;  // IAF -> wigner transform (complex)

    void wigner_distribution(const SimState& state);

    void husimi_distribution(const SimState& state);

   public:
    PhaseSpaceDistribution(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(PhaseSpaceDistribution)
};

class Potential : public ObservableFunctor {
   public:
    Potential(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(Potential)
};

class WaveFunction : public ObservableFunctor {
    double t_prev;
    blaze::DynamicMatrix<double> psi_re_im;

   public:
    WaveFunction(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(WaveFunction)
};

class Energy : public ObservableFunctor {
    const Cosmology& cosmo;
    const int N;
    const double L;
    const double dx;
    blaze::CompressedMatrix<double> grad;
    blaze::DynamicVector<double> energies;

   public:
    Energy(const Parameters& p, const Cosmology& cosmo_);
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(Energy)
};

class Entropy : public ObservableFunctor {
    const Parameters p;
    const Cosmology cosmo;
    const int N;
    const double dx;
    const double dk;
    double t_prev;
    blaze::DynamicVector<double> S;
    std::unique_ptr<ObservableFunctor> phasespace;

   public:
    Entropy(const Parameters& p, const Cosmology&);
    ObservableFunctor::ReturnType compute(
        const SimState& state,
        const std::unordered_map<
            std::string, std::unique_ptr<ObservableFunctor>>& obs) override;
    REGISTER(Entropy)
};

}  // namespace Observable
#endif
