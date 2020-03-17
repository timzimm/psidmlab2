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
#include "interaction/poisson_fft.h"
#include "parameters.h"
#include "parameters_fwd.h"
#include "state.h"

using namespace blaze;

namespace {
// This test suite tests the implementation of the Poisson FFT class with
// finite range
//
//              del^2_x V(x) = f(x) + PBC
//
class PoissonFFTFixture
    : public testing::TestWithParam<std::tuple<int, double>> {
   protected:
    PoissonFFTFixture()
        // Part of parameter file required to carry out the test
        : p({{"Domain",
              {{"N", std::get<0>(GetParam())}, {"L", std::get<1>(GetParam())}}},
             {"Simulation", {{"interaction", {{"epsilon", 1}}}}}}),
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
    Poisson::FFT solver;
    blaze::DynamicVector<double> f;
    blaze::DynamicVector<double> V_ref;
};

TEST_P(PoissonFFTFixture, AnalyticalvsInPlace) {
    const double tolerance = 1e-6;
    solver.solve(f, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(f - V_ref), 0, tolerance);
}

TEST_P(PoissonFFTFixture, AnalyticalvsOutOfPlace) {
    const double tolerance = 1e-6;
    DynamicVector<double> V(box.N);
    DynamicVector<double> f_copy = f;
    solver.solve(V, f);
    // Result correct?
    EXPECT_NEAR(maxNorm(V - V_ref), 0, tolerance);
    // Source altered?
    EXPECT_EQ(f_copy, f);
}

INSTANTIATE_TEST_SUITE_P(PowerOf2N, PoissonFFTFixture,
                         testing::Combine(testing::Values(256, 512, 1024, 2048,
                                                          4096),
                                          testing::Values(1.0, 2 * M_PI, 13)));

INSTANTIATE_TEST_SUITE_P(OddN, PoissonFFTFixture,
                         testing::Combine(testing::Values(15, 17, 193, 1881),
                                          testing::Values(1.0, 2 * M_PI, 13)));

}  // namespace
