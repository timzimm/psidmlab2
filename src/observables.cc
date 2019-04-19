#include "observables.h"
#include "blaze/math/UniformVector.h"
#include "parameters.h"
#include "state.h"

namespace Observable {

// TODO The initializer list is too messy
DensityContrast::DensityContrast(const Parameters& p)
    : sigma_x(p["Analysis"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      ws(husimi ? p["Analysis"]["linear_convolution"].get<bool>() : 0,
         husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(N_kernel) {
    // Kernel construction in x-space
    for (int i = 0; i < N_kernel; ++i) {
        gaussian_kernel[i] = (i - N_kernel / 2) * dx;
    }
    gaussian_kernel =
        exp(-0.5 / (sigma_x * sigma_x) * gaussian_kernel * gaussian_kernel);
    gaussian_kernel /= blaze::sum(gaussian_kernel);
}

blaze::DynamicMatrix<double, blaze::columnMajor> DensityContrast::compute(
    const SimState& state) {
    auto psi2 = blaze::real(state.psis % state.psis);
    blaze::UniformVector<double> one(psi2.rows());
    blaze::DynamicVector<double> delta = psi2 * state.lambda - one;

    if (husimi) {
        discrete_convolution(ws, gaussian_kernel, delta);
        delta = blaze::subvector(ws.signal_padded, ws.N_offset, N);
    }

    return blaze::expand(delta, 1);
}

}  // namespace Observable
