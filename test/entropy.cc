#include <gtest/gtest.h>
#include <boost/variant.hpp>
#include <fstream>
#include "gtest/gtest-param-test.h"

#include "blaze/math/expressions/DVecNormExpr.h"
#include "blaze/math/functors/L2Norm.h"
#include "convolution_functions.h"
#include "cosmology.h"
#include "interfaces.h"
#include "parameters.h"
#include "parameters_fwd.h"
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
// and computed via the convolution theorem.
class EntropyTest
    : public testing::TestWithParam<std::tuple<int, double, double>> {
   protected:
    void SetUp() override {
        std::tie(N, L, sigma_x) = GetParam();
        dx = L / N;
        // Setup part of json that is required to init state and observable
        std::ofstream a_of_t("a_of_t.txt");
        a_of_t << "0\t0\n";
        a_of_t.close();

        p["Simulation"]["N"] = N;
        p["Simulation"]["L"] = L;

        p["Cosmology"]["model"] = 1;
        p["Cosmology"]["scalefactor_file"] = "./a_of_t.txt";

        p["Observables"]["PhaseSpaceDistribution"]["linear_convolution"] =
            linear;
        p["Observables"]["PhaseSpaceDistribution"]["sigma_x"] = sigma_x;
        p["Observables"]["PhaseSpaceDistribution"]["patch"] = {};
        p["Observables"]["Entropy"]["compute_at"] = {0};

        cosmo = Cosmology(p["Cosmology"]);
        obs["Entropy"] =
            ObservableFunctor::make("Observable::Entropy", p, cosmo);
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

// It can be shown that
//                                 S(f_H) >= 1
// for any f_H. We test this by computing S for a couple of smoothed random psis
// which are localized and have unit norm.
TEST_P(EntropyTest, WehrlsConjecture) {
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
    ASSERT_NEAR(dx * std::pow(real(l2Norm(state.psi)), 2), 1, 1e-14);

    ObservableFunctor* entropy = obs["Entropy"].get();
    auto result = entropy->compute(state, obs);
    auto S = boost::get<const DynamicVector<double>&>(result);
    EXPECT_GE(S[0], 1);
}

// For psi(x) as defined below we have:
//                  S = - int dx dk f log f = 2 pi
// analytically.
TEST_P(EntropyTest, GaussPsiEntropyIsTwoPi) {
    state.psi.resize(N);
    auto x = linspace(N, -L / 2, L / 2 - L / N);
    state.psi = std::pow(2 * M_PI * sigma_x * sigma_x, -0.25) *
                exp(-x * x / (4 * sigma_x * sigma_x));

    ObservableFunctor* entropy = obs["Entropy"].get();
    auto result = entropy->compute(state, obs);
    auto S = boost::get<const DynamicVector<double>&>(result);
    EXPECT_NEAR(S[0], 2 * M_PI, 1e-4);
}

INSTANTIATE_TEST_SUITE_P(
    Wehrl, EntropyTest,
    testing::Combine(testing::Values(1024, 1025, 2048, 2049),
                     testing::Values(100, 200), testing::Values(0.25, 0.5, 1)));
}  // namespace
