#include <gtest/gtest.h>
#include <boost/variant.hpp>
#include <fstream>
#include <limits>

#include "convolution_functions.h"  //for psi construction
#include "cosmology.h"
#include "domain.h"
#include "interfaces.h"  //we want to test observables
#include "parameters.h"
#include "state.h"

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
    HusimiTest()
        // Part of parameter file required to carry out the test
        : sigma_x(std::get<2>(GetParam())),
          p({{"Simulation",  // Remove when everything uses Domain
              {
                  {"N", std::get<0>(GetParam())},
                  {"L", std::get<1>(GetParam())},
              }},
             {"Domain",
              {
                  {"N", std::get<0>(GetParam())},
                  {"L", std::get<1>(GetParam())},
              }},
             {"Cosmology",
              {
                  {"model", 1},
                  {"scalefactor_file", "./a_of_t.txt"},
              }},
             {"Observables",
              {{"PhaseSpaceDistribution",
                {{"linear_convolution", false},
                 {"sigma_x", sigma_x},
                 {"compute_at", {0}},
                 {"patch", {}}}}}}}),
          cosmo(p),
          box(p) {
        // Write a_of_t file to directory
        std::ofstream a_of_t("a_of_t.txt");
        a_of_t << "0\t0\n";
        a_of_t.close();

        // Allocate the entropy observable
        obs["PhaseSpaceDistribution"] = ObservableFunctor::make(
            "Observable::PhaseSpaceDistribution", p, cosmo);
    }

    const double sigma_x;
    const Parameters p;
    const Cosmology cosmo;
    const Domain box;

    SimState state;
    ObservableMap obs;
};

// Husimis distribution is non-negative definite and bounded. For unit norm psi
// one can show:
//                              0 <= f <= 1/pi
// for any f_H. Since we add a small offset to f this becomes:
//                            eps <= f <= 1/pi + eps
// We test this by computing f for a couple of smoothed random psis
// which are localized and have unit norm.
TEST_P(HusimiTest, IsBounded) {
    const double eps = std::numeric_limits<double>::min();
    const int N_kernel = 10;
    const int N = box.N;
    const double dx = box.dx;

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

    auto x = linspace(box.N, box.xmin, box.xmax);
    auto k = linspace(box.N, box.kmin, box.kmax);

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

// For our application psi is noramlized as follows
//              int dx |psi^2| = L
// Check if the Husimi function obeys the same normalization
// if integrated over x and k, which it has to analytically.
TEST_P(HusimiTest, NormalizationAsPsi) {
    const double tolerance = 1e-6;

    auto x = linspace(box.N, box.xmin, box.xmax);
    auto k = linspace(box.N, box.kmin, box.kmax);

    state.psi = std::pow(2 * M_PI * sigma_x * sigma_x, -0.25) *
                exp(-x * x / (4 * sigma_x * sigma_x));
    double norm =
        sqrt(box.L) / sqrt(box.dx * sum(real(conj(state.psi) * state.psi)));
    state.psi *= norm;

    // int dx |psi^2| = L ?
    EXPECT_NEAR(box.L, box.dx * sum(real(conj(state.psi) * state.psi)),
                tolerance);

    ObservableFunctor* husimi = obs["PhaseSpaceDistribution"].get();
    auto result = husimi->compute(state, obs);
    const DynamicMatrix<double>& f =
        boost::get<const DynamicMatrix<double>&>(result);
    double norm_f = 1.0 / (2 * M_PI) * box.dx * box.dk * sum(f);

    // int dx dk/2pi fH = L ?
    EXPECT_NEAR(norm_f, box.L, tolerance);
}

INSTANTIATE_TEST_SUITE_P(PhaseSpace, HusimiTest,
                         testing::Combine(testing::Values(1024, 1025, 2048,
                                                          2049),
                                          testing::Values(100, 200),
                                          testing::Values(0.25, 0.5, 0.75)));
}  // namespace
