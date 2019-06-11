#include "observables.h"
#include "blaze/math/Columns.h"
#include "blaze/math/Elements.h"
#include "blaze/math/Submatrix.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>

namespace Observable {

DensityContrast::DensityContrast(const Parameters& p)
    : sigma_x(p["Observables"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Observables"]["linear_convolution"].get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      t_prev(-1),
      ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(husimi ? N_kernel : 0),
      delta(N) {
    if (husimi) {
        // Kernel construction in x-space
        for (int i = 0; i < N_kernel; ++i) {
            gaussian_kernel[i] = (i - N_kernel / 2) * dx;
        }
        gaussian_kernel =
            exp(-0.5 / (sigma_x * sigma_x) * gaussian_kernel * gaussian_kernel);
        gaussian_kernel /= blaze::sum(gaussian_kernel);
    }
}

ObservableFunctor::ReturnType DensityContrast::compute(const SimState& state) {
    using namespace blaze;
    if (t_prev < state.tau) {
        delta = delta_from(state);

        if (husimi) {
            discrete_convolution(ws, gaussian_kernel, delta);
            delta = blaze::subvector(ws.signal_padded, ws.N_offset, N);
            if (!linear) {
                // TODO Benchmark against element selection
                std::rotate(std::begin(delta), std::begin(delta) + N_kernel / 2,
                            std::end(delta));
            }
        }
        t_prev = state.tau;
    }

    return delta;
}

PhaseSpaceDistribution::PhaseSpaceDistribution(const Parameters& p)
    : sigma_x(std::sqrt(2) * p["Observables"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Observables"]["linear_convolution"].get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      L(p["Simulation"]["L"].get<double>()),
      dx(L / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      t_prev(-1),
      ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
      gaussian_kernel(0),
      f(N, N, 0),
      iaf(N, N, 0),
      c2c(nullptr) {
    if (husimi) {
        gaussian_kernel.resize(N_kernel);
        // Kernel construction in x-space
        for (int i = 0; i < N_kernel; ++i) {
            gaussian_kernel[i] = (i - N_kernel / 2) * dx;
        }
        gaussian_kernel =
            exp(-0.5 / (sigma_x * sigma_x) * gaussian_kernel * gaussian_kernel);
        gaussian_kernel /= blaze::sum(gaussian_kernel);
    }
    auto iaf_ptr = reinterpret_cast<fftw_complex*>(iaf.data());
    c2c = fftw_plan_many_dft(1, &N, N, iaf_ptr, nullptr, 1, iaf.spacing(),
                             iaf_ptr, nullptr, 1, iaf.spacing(), FFTW_FORWARD,
                             FFTW_ESTIMATE);
}

void PhaseSpaceDistribution::wigner_distribution(const SimState& state) {
    using namespace blaze;
    for (int m = 0; m < state.M; ++m) {
        auto psi = column(state.psis, m);
        for (int i = 0; i < N; ++i) {
            // TODO Do we need to reinit to zero on each new iteration?
            auto iaf_i = column(iaf, i);
            int min = std::min(i, N - 1 - i);
            int total = 2 * min + 1;
            auto lag_plus =
                elements(psi, [&](int tau) { return i + (tau - min); }, total);
            auto lag_minus =
                elements(psi, [&](int tau) { return i - (tau - min); }, total);
            auto index = elements(
                iaf_i, [&](int tau) { return (N + tau - min) % N; }, total);
            index = lag_plus * conj(lag_minus);
            /* for (int j = 0; 2 * j < N; ++j) { */
            /*     iaf_i[2 * j + 1] *= -1; */
            /* } */
        }
        fftw_execute(c2c);

        f += state.lambda[m] * real(iaf);
    }
}

void PhaseSpaceDistribution::husimi_distribution(const SimState& state) {
    using namespace blaze;

    DynamicMatrix<std::complex<double>, columnMajor> husimi_fft(N, N);

    double normal = L / (N * sqrt(4 * M_PI * sigma_x * sigma_x));

    // Build up blaze expression representing the gaussian kernel
    DynamicVector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    auto idx_mat = expand(idx, N);
    auto diff = trans(idx_mat) - idx_mat;
    auto gauss =
        declsym(normal * exp(-1.0 / (4 * N_kernel * N_kernel) * diff % diff));

    // DFT to obtain psi_husimi
    fftw_complex* in_out = reinterpret_cast<fftw_complex*>(husimi_fft.data());
    auto plan = fftw_plan_many_dft(
        1, &N, N, in_out, nullptr, 1, husimi_fft.spacing(), in_out, nullptr, 1,
        husimi_fft.spacing(), FFTW_FORWARD, FFTW_ESTIMATE);

    // Construct mixed state distribution pure state by pure state
    for (int m = 0; m < state.M; ++m) {
        auto psi = column(state.psis, m);
        husimi_fft = gauss % expand(psi, N);
        fftw_execute(plan);
        f += state.lambda[m] * real(husimi_fft % husimi_fft);
    }

    fftw_destroy_plan(plan);

    // Move zero velocity into matrix center.
    for (int i = 0; i < N; ++i) {
        for (int j = 0; 2 * j < N; ++j) {
            f(i, 2 * j + 1) *= -1;
        }
    }
}

ObservableFunctor::ReturnType PhaseSpaceDistribution::compute(
    const SimState& state) {
    if (t_prev < state.tau) {
        husimi ? husimi_distribution(state) : wigner_distribution(state);
        t_prev = state.tau;
    }
    return f;
}

PhaseSpaceDistribution::~PhaseSpaceDistribution() { fftw_destroy_plan(c2c); }

Potential::Potential(const Parameters& p){};

inline ObservableFunctor::ReturnType Potential::compute(const SimState& state) {
    return state.V;
}

WaveFunction::WaveFunction(const Parameters& p){};

inline ObservableFunctor::ReturnType WaveFunction::compute(
    const SimState& state) {
    return state.psis;
}

}  // namespace Observable

