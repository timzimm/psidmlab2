#include <gtest/gtest.h>
#include <boost/variant.hpp>
#include <fstream>
#include "gtest/gtest-param-test.h"

#include "convolution_functions.h"  //for psi construction
#include "cosmology.h"
#include "domain.h"
#include "interfaces.h"
#include "parameters.h"
#include "state.h"

using namespace blaze;

using ObservableMap =
    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>;

namespace {
// This test suite tests the implementation of Wehrl's entropy defined as:
//
//                      S[f] = - int dx dk f log f
//
// with f as Husimi function.
class EntropyTest
    : public testing::TestWithParam<std::tuple<int, double, double>> {
   protected:
    EntropyTest()
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
        obs["Entropy"] =
            ObservableFunctor::make("Observable::Entropy", p, cosmo);
    }
    const double sigma_x;
    const Parameters p;
    const Cosmology cosmo;
    const Domain box;

    SimState state;
    ObservableMap obs;
};

// It can be shown that
//                                 S(f_H) >= 1
// for any f_H. We test this by computing S for a couple of smoothed random psis
// which are localized and have unit norm.
TEST_P(EntropyTest, WehrlsConjecture) {
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

    ObservableFunctor* entropy = obs["Entropy"].get();
    auto result = entropy->compute(state, obs);
    const DynamicVector<double>& S =
        boost::get<const DynamicVector<double>&>(result);

    EXPECT_GE(S[0], 1);
}

// For psi(x) as defined below we have:
//                  S = - int dx dk f log f = 2 pi
// analytically.
TEST_P(EntropyTest, GaussPsiEntropyIsTwoPi) {
    const double tolerance = 1e-6;
    const int N = box.N;
    const double L = box.L;

    state.psi.resize(N);
    auto x = linspace(N, -L / 2, L / 2 - L / N);
    state.psi = std::pow(2 * M_PI * sigma_x * sigma_x, -0.25) *
                exp(-x * x / (4 * sigma_x * sigma_x));

    ObservableFunctor* entropy = obs["Entropy"].get();
    auto result = entropy->compute(state, obs);
    const DynamicVector<double>& S =
        boost::get<const DynamicVector<double>&>(result);

    EXPECT_NEAR(S[0], 2 * M_PI, tolerance);
}

INSTANTIATE_TEST_SUITE_P(Wehrl, EntropyTest,
                         testing::Combine(testing::Values(1024, 1025, 2048,
                                                          2049),
                                          testing::Values(100, 200),
                                          testing::Values(0.25, 0.5, 0.75)));
}  // namespace
