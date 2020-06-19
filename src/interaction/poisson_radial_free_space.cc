#include "interaction/poisson_radial_free_space.h"
#include "fftw.h"
#include "parameters.h"
#include "state.h"

using namespace blaze;

namespace Poisson {
RadialFreeSpace::RadialFreeSpace(const Parameters &p, const SimState &state)
    : box(p), dst(nullptr), source_k(box.N) {
    if (box.bc != Domain::BoundaryCondition::HomogeneousDirichlet) {
        std::cout << ERRORTAG("Wrong Boundary Conditions") << std::endl;
        exit(1);
    }
    // FFTW only stores non-redundant modes. Using
    //              phi_i  = phi_N+i    (periodicity) (i)
    //              phi_i  = -phi_(N-i) (oddness)     (ii)
    // This implies only the first box.N integers with DFT_size=2(box.N+1) are
    // non-redundant. box.N starts counting from i=1 because i=0 must be zero
    // due to both conditions above. We use the DST to enforce homogeneous
    // Dirichlet at i=0 and i=box.N+1.
    // ________________________ DFT_size _______________________
    // |                                                        |       x0 (i)
    // x0       x1      x2      x3      x4      x5      x6      x7      x8
    // 0       |                  |    -x4     -x3     -x2     -x1     -x0 (ii)
    //         ------ box.N -------
    //
    // Transformations are most efficient if DFT_size is a power of two
    // implying box.N = 2*k - 1. However ANY box.N can be used
    //
    // Also note that RODFT00 is its own inverse. Hence one plan is enough.
    dst = make_fftw_plan_r2r_1d(box.N, source_k.data(), source_k.data(),
                                FFTW_RODFT00, FFTW_MEASURE);
}

void RadialFreeSpace::solve(SimState &state) {
    solve(state.V, rho_from(state));
}

void RadialFreeSpace::solve(blaze::DynamicVector<double> &V,
                            const blaze::DynamicVector<double> &source) {
    const int M = box.N + 1;
    static auto r = linspace(box.N, box.xmin, box.xmax);

    // DST is unnormalized
    source_k = 1.0 / (M * r) * source;
    fftw_execute(dst.get());
    auto kx = linspace(box.N, box.kmin, box.kmax);
    source_k /= -1.0 * kx * kx;
    fftw_execute(dst.get());
    V = 1.0 / (4 * M_PI) * (1.0 / r * source_k + 1.0 / box.L);
}

}  // namespace Poisson
