#include <blaze/math/dense/DynamicVector.h>
#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>
#include <fstream>
#include <iterator>
#include <limits>
#include <vector>

#include "domain.h"
#include "interaction/poisson_fft_epsilon.h"
#include "parameters.h"
#include "parameters_fwd.h"
#include "state.h"

using namespace blaze;

using ObservableMap =
    std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>;

namespace {
// This test suite tests the implementation of the Poisson FFT class with
// finite range
//
//              del^2_x V(x) - epsilon^2 V(x) = f(x) + PBC
//
class PoissonFFTEpsilonFixture
    : public testing::TestWithParam<std::tuple<int, double>> {
   protected:
    PoissonFFTEpsilonFixture()
        // Part of parameter file required to carry out the test
        : p({{"Domain",
              {{"N", std::get<0>(GetParam())}, {"L", std::get<1>(GetParam())}}},
             {"Simulation", {{"interaction", {{"epsilon", 1}}}}}}),
          box(p),
          solver(p, state),
          f(box.N, 0),
          V_ref(box.N, 0) {
        // For epsilon = 1, x \in [0,L) and:
        // f(x) = -exp^(sin(2pi/L * x)*(4pi^2/L^2*sin(x)-4pi^2/L^2*cos^2(x)+1)
        // we have:
        //          V(x) = exp^(sin(2pi/L * x))
        // analytically.
        double omega = 2 * M_PI / box.L;
        double omega2 = omega * omega;
        auto sinx = sin(omega * linspace(box.N, box.xmin, box.xmax));
        auto cosx = cos(omega * linspace(box.N, box.xmin, box.xmax));
        f = -exp(sinx) * (omega2 * sinx - omega2 * cosx * cosx + 1);
        V_ref = exp(sinx);
    }

    const Parameters p;
    const Domain box;

    SimState state;
    Poisson::FFTEpsilon solver;
    blaze::DynamicVector<double> f;
    blaze::DynamicVector<double> V_ref;
};  // namespace

// Check if error decays as N increases and write the error to file.
TEST(PoissonFFTEpsilonTest, ErrorDecays) {
    std::ofstream error_file = std::ofstream("errors.txt");
    std::vector<double> errors;
    for (int i = 1; i < 12; ++i) {
        const int N = 1 << i;
        Parameters p = {{"Domain", {{"N", N}, {"L", 2 * M_PI}}},
                        {"Simulation", {{"interaction", {{"epsilon", 1}}}}}};
        SimState state;
        Domain box(p);
        Poisson::FFTEpsilon solver(p, state);
        auto sinx = sin(linspace(box.N, box.xmin, box.xmax));
        DynamicVector<double> f = -exp(sinx) * (sinx * sinx + sinx);
        DynamicVector<double> V_ref = exp(sinx);
        solver.solve(f, f);
        double error = maxNorm(f - V_ref);
        // error decaying compared to last run with smaller N?
        if (error > 10 * std::numeric_limits<double>::epsilon() &&
            errors.size() > 0)
            EXPECT_LT(error, errors.back());

        errors.push_back(error);
        error_file << box.N << "\t" << error << std::endl;
    }
    error_file.close();
}

TEST_P(PoissonFFTEpsilonFixture, AnalyticalvsInPlace) {
    const double tolerance = 1e-6;
    solver.solve(f, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(f - V_ref), 0, tolerance);
}

TEST_P(PoissonFFTEpsilonFixture, AnalyticalvsOutOfPlace) {
    const double tolerance = 1e-6;
    DynamicVector<double> V(box.N);
    DynamicVector<double> f_copy = f;
    solver.solve(V, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(V - V_ref), 0, tolerance);
    // Source altered?
    EXPECT_EQ(f_copy, f);
}

INSTANTIATE_TEST_SUITE_P(PowerOf2N, PoissonFFTEpsilonFixture,
                         testing::Combine(testing::Values(256, 512, 1024, 2048,
                                                          4096),
                                          testing::Values(1.0, 2 * M_PI, 13)));

INSTANTIATE_TEST_SUITE_P(OddN, PoissonFFTEpsilonFixture,
                         testing::Combine(testing::Values(15, 17, 193, 1881),
                                          testing::Values(1.0, 2 * M_PI, 13)));

}  // namespace
