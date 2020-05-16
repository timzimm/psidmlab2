#include <blaze/math/dense/DynamicVector.h>
#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "blaze/math/functors/Max.h"
#include "domain.h"
#include "interaction/periodic_convolution.h"
#include "parameters.h"
#include "parameters_fwd.h"
#include "state.h"

using namespace blaze;

namespace {
// This test suite tests the implementation of the Periodic Convolution Class
class PoissonFixture : public testing::TestWithParam<std::tuple<int, double>> {
   protected:
    PoissonFixture()
        // Part of parameter file required to carry out the test
        : p({{"Domain",
              {{"N", std::get<0>(GetParam())}, {"L", std::get<1>(GetParam())}}},
             {"Simulation", {{"interaction", {{"name", "Poisson"}}}}}}),
          box(p),
          solver(p, state),
          f(box.N, 0),
          V_ref(box.N, 0) {
        // For epsilon = 1, x \in [0,L) and:
        //      f(x) = 4*pi^2/L^2 sin(2pi/L * x)
        // we have:
        //          V(x) = -sin(2pi/L * x)
        // analytically.
        double omega = 2 * M_PI / box.L;
        auto sinx = sin(omega * linspace(box.N, box.xmin, box.xmax));
        f = omega * omega * sinx;
        V_ref = -sinx;
    }

    const Parameters p;
    const Domain box;

    SimState state;
    PeriodicConvolution solver;
    blaze::DynamicVector<double> f;
    blaze::DynamicVector<double> V_ref;
};

TEST_P(PoissonFixture, AnalyticalvsInPlace) {
    const double tolerance = 1e-6;
    solver.solve(f, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(f - V_ref), 0, tolerance);
}

TEST_P(PoissonFixture, AnalyticalvsOutOfPlace) {
    const double tolerance = 1e-6;
    DynamicVector<double> V(box.N);
    DynamicVector<double> f_copy = f;
    solver.solve(V, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(V - V_ref), 0, tolerance);
    // Source altered?
    EXPECT_EQ(f_copy, f);
}

class ScreenedPoissonFixture
    : public testing::TestWithParam<std::tuple<int, double>> {
   protected:
    ScreenedPoissonFixture()
        // Part of parameter file required to carry out the test
        : p({{"Domain",
              {{"N", std::get<0>(GetParam())}, {"L", std::get<1>(GetParam())}}},
             {"Simulation",
              {{"interaction",
                {{"name", "ScreenedPoisson"}, {"epsilon", 1}}}}}}),
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
    PeriodicConvolution solver;
    blaze::DynamicVector<double> f;
    blaze::DynamicVector<double> V_ref;
};  // namespace

TEST_P(ScreenedPoissonFixture, AnalyticalvsInPlace) {
    const double tolerance = 1e-6;
    solver.solve(f, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(f - V_ref), 0, tolerance);
}

TEST_P(ScreenedPoissonFixture, AnalyticalvsOutOfPlace) {
    const double tolerance = 1e-6;
    DynamicVector<double> V(box.N);
    DynamicVector<double> f_copy = f;
    solver.solve(V, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(V - V_ref), 0, tolerance);
    // Source altered?
    EXPECT_EQ(f_copy, f);
}

INSTANTIATE_TEST_SUITE_P(PowerOf2N, PoissonFixture,
                         testing::Combine(testing::Values(256, 512, 1024, 2048,
                                                          4096),
                                          testing::Values(1.0, 2 * M_PI, 13)));

INSTANTIATE_TEST_SUITE_P(OddN, PoissonFixture,
                         testing::Combine(testing::Values(15, 17, 193, 1881),
                                          testing::Values(1.0, 2 * M_PI, 13)));

INSTANTIATE_TEST_SUITE_P(PowerOf2N, ScreenedPoissonFixture,
                         testing::Combine(testing::Values(256, 512, 1024, 2048,
                                                          4096),
                                          testing::Values(1.0, 2 * M_PI, 13)));

INSTANTIATE_TEST_SUITE_P(OddN, ScreenedPoissonFixture,
                         testing::Combine(testing::Values(15, 17, 193, 1881),
                                          testing::Values(1.0, 2 * M_PI, 13)));

}  // namespace
