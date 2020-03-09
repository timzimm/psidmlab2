#include <gtest/gtest.h>
#include <boost/variant.hpp>
#include <fstream>
#include <limits>

#include "blaze/math/dense/DynamicMatrix.h"
#include "convolution_functions.h"
#include "cosmology.h"
#include "observables_common.h"

using namespace blaze;

using ObservableMap =
    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>;

namespace {
// This test suite tests the implementation of the Husimi Distribution defined
// as:
//
// f(x,k) =
//   | 1/(2*pi*s^2)^(1/4) int dx' dk exp(-(x-x')^2/(4s^2) - ikx') psi(x') |^2
//
class HusimiTest
    : public testing::TestWithParam<std::tuple<int, double, double>> {
   protected:
    void SetUp() override {
        std::tie(N, L, sigma_x) = GetParam();
        dx = L / N;

        // Box setup
        p["Simulation"]["N"] = N;
        p["Simulation"]["L"] = L;
        state.psi.resize(N);

        // Cosmology setup
        p["Cosmology"]["model"] = 1;
        p["Cosmology"]["scalefactor_file"] = "./a_of_t.txt";
        std::ofstream a_of_t("a_of_t.txt");
        a_of_t << "0\t0\n";
        a_of_t.close();
        cosmo = Cosmology(p["Cosmology"]);

        // observable setup
        const std::string name = "PhaseSpaceDistribution";
        p["Observables"][name]["compute_at"] = {0};
        p["Observables"][name]["linear_convolution"] = linear;
        p["Observables"][name]["sigma_x"] = sigma_x;
        p["Observables"][name]["patch"] = {};
        obs[name] = ObservableFunctor::make("Observable::" + name, p, cosmo);
    }
    Parameters p;
    SimState state;
    Cosmology cosmo;
    ObservableMap obs;
    int N;
    double L;
    double dx;
    double sigma_x;
    bool linear;
};

// Husimis distribution is non-negative definite and bounded. For unit norm psi
// one can show:
//                              0 <= f <= 1/pi
// for any f_H. Since we add a small offset to f this becomes:
//                            eps <= f <= 1/pi + eps
// We test this by computing f for a couple of smoothed random psis
// which are localized and have unit norm.
TEST_P(HusimiTest, fIsBounded) {
    const double eps = std::numeric_limits<double>::min();
    const int N_kernel = 10;
    convolution_ws<std::complex<double>> ws(false, N, N_kernel);

    // Simple box filter of length N_kernel * dx
    ws.kernel_padded = 1;

    // Randomize psi
    auto psi_center = subvector(ws.signal_padded, N / 4, N / 2);
    for (auto& psi : psi_center) {
        psi = blaze::rand<double>(0.0, 1.0);
    }

    // Smooth psi
    discrete_convolution(ws);

    // Establish unit L2-norm
    double norm = std::sqrt(dx) * std::real(l2Norm(ws.signal_padded));
    state.psi = ws.signal_padded / norm;
    ASSERT_NEAR(std::sqrt(dx) * std::real(l2Norm(state.psi)), 1, 1e-14);

    ObservableFunctor* husimi = obs["PhaseSpaceDistribution"].get();
    auto result = husimi->compute(state, obs);
    const DynamicMatrix<double>& f =
        boost::get<const DynamicMatrix<double>&>(result);

    EXPECT_LE(eps, min(f));
    EXPECT_LE(max(f), eps + 1 / M_PI);
}

// For psi(x) as defined below we have:
//          f(x,k) = exp(-(sigma_x*k)^2*exp(-x^2/(4sigma^2))
// analytically.
TEST_P(HusimiTest, IsGaussian) {
    const double tolerance = 1e-6;
    const double dx = L / N;
    const double dk = 2 * M_PI / L;

    // Total x-k grid (see dicussion in observables/phasespace_distribution.cc)
    const double xmin = -L / 2;
    const double xmax = -xmin - dx;
    const double kmin = -dk * (N / 2);
    const double kmax = (N % 2) ? -kmin : -kmin - dk;

    auto x = linspace(N, xmin, xmax);
    auto k = linspace(N, kmin, kmax);

    state.psi = std::pow(2 * M_PI * sigma_x * sigma_x, -0.25) *
                exp(-x * x / (4 * sigma_x * sigma_x));

    // Reference solution (see above)
    auto f_ref = exp(-sigma_x * sigma_x * k * k) *
                 trans(exp(-0.25 / (sigma_x * sigma_x) * x * x));

    ObservableFunctor* husimi = obs["PhaseSpaceDistribution"].get();
    auto result = husimi->compute(state, obs);
    const DynamicMatrix<double>& f =
        boost::get<const DynamicMatrix<double>&>(result);
    EXPECT_NEAR(maxNorm(f - f_ref), 0, tolerance);
}

INSTANTIATE_TEST_SUITE_P(PhaseSpace, HusimiTest,
                         testing::Combine(testing::Values(1024, 1025, 2048,
                                                          2049),
                                          testing::Values(100, 200),
                                          testing::Values(0.25, 0.5, 0.75)));
}  // namespace
