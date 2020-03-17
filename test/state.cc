#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>
#include <fstream>
#include <limits>
#include <vector>

#include "convolution_functions.h"
#include "domain.h"
#include "parameters.h"
#include "state.h"

using namespace blaze;

namespace {
// This test suite tests the implementation of the Simulationn State struct
//
class SimStateFixture : public testing::TestWithParam<std::tuple<int, double>> {
   protected:
    SimStateFixture()
        // Part of parameter file required to carry out the test
        : p({
              {"Domain",
               {{"N", std::get<0>(GetParam())},
                {"L", std::get<1>(GetParam())}}},
          }),
          box(p) {
        state.psi.resize(box.N);
        state.V.resize(box.N);
    }

    const Parameters p;
    const Domain box;

    SimState state;
};

// Parselvals Theorem
TEST_P(SimStateFixture, Parseval) {
    const double tolerance = 1e-10;
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

    // Establish L2-norm
    double norm = std::sqrt(dx) * std::real(l2Norm(ws.signal_padded));
    state.psi = sqrt(box.L) * ws.signal_padded / norm;

    ASSERT_NEAR(dx * std::real(sum(conj(state.psi) * state.psi)), box.L,
                tolerance);
    double norm2_real = N;

    for (int i = 0; i < 100000; ++i) {
        state.transform(SimState::Representation::Position);
        state.transform(SimState::Representation::Momentum);
        double norm2_fft = sum(real(conj(state.psi) * state.psi));
        EXPECT_NEAR(N * norm2_fft, norm2_real, tolerance);
    }
}

INSTANTIATE_TEST_SUITE_P(StateTest, SimStateFixture,
                         testing::Combine(testing::Values(256, 255, 345, 512),
                                          testing::Values(1.0, 4.0, 100.0)));

}  // namespace
