#include "interaction/screened_poisson_fc.h"
#include "fftw.h"
#include "fftw3.h"
#include "parameters.h"
#include "state.h"

using namespace blaze;

namespace ScreenedPoisson {
FastConvolution::FastConvolution(const Parameters &p, const SimState &state)
    : box(p),
      epsilon(p["Simulation"]["interaction"]["epsilon"]),
      dst(nullptr),
      source_k(box.N) {
    // FFTW only stores non-redundant modes. Using
    //              phi_i  = phi_N+i    (periodicity)
    //              phi_i  = -phi_(N-i) (oddness)
    // This implies only the first n integers with N=2(n+1) are non-redundant.
    // n starts counting from i=1 because i=0 must be zero due to both
    // conditions above. We use the DST to enforce homogeneous Dirichlet at
    // i=0 and i=N+1. To satisfy the second root we must have N, i.e. number
    // of non-redundant elements on the grid, as odd:
    // ________________________ DFT size _______________________
    // |                                                        |
    // x0       x1      x2      x3      x4      x5      x6      x7
    // 0       |                  |    -x4     -x3     -x2     -x1
    //         -------  N  --------     ^
    //
    if (box.N % 2 == 0) {
        std::cerr
            << ERRORTAG(
                   "N must be odd to enforce homogenous Dirchlet condition")
            << std::endl;
        exit(1);
    }
    // Also note that RODFT00 is its own inverse. Hence one plan is enough.
    dst = make_fftw_plan_r2r_1d(box.N, source_k.data(), source_k.data(),
                                FFTW_RODFT00, FFTW_ESTIMATE);
}

void FastConvolution::solve(SimState &state) {
    solve(state.V, rho_from(state));
}

// Fast Convolution Algorithm as introduced in:
// Y.Zhang, X.Dong Journal of Compuational Physics 230 (2011) 2660-2676
void FastConvolution::solve(blaze::DynamicVector<double> &V,
                            const blaze::DynamicVector<double> &source) {
    const int M = box.N + 1;

    // DST is unnormalized
    source_k = 1.0 / M * source;
    fftw_execute(dst.get());

    // sum index
    auto l = linspace(M - 1, 1, M - 1);
    // wavenumber
    auto mu_l = M_PI * l / box.L;
    // x coordinates of non-boundary points
    auto x_j = box.xmin + box.L / M * l;
    // free space Greens function
    auto G_l = -1.0 / (epsilon * epsilon + mu_l * mu_l);

    // Uniform sums, i.e. independent of x can be computed once
    double S1 = sum(mu_l / (2 * epsilon) * G_l * source_k);
    double S2 = sum(mu_l / (2 * epsilon) * G_l * cos(M_PI * l) * source_k);

    // x dependent sum is accelerated with DST
    source_k *= 0.5 * G_l;
    fftw_execute(dst.get());

    V = S1 * exp(epsilon * (box.xmin - x_j)) -
        S2 * exp(epsilon * (x_j - box.xmax)) + source_k;
}

}  // namespace ScreenedPoisson
