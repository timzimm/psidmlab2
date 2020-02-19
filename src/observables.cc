#include "observables.h"
#include "cosmology.h"
#include "parameters.h"
#include "state.h"

#include <blaze/math/Columns.h>
#include <blaze/math/Elements.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Submatrix.h>
#include <algorithm>

using namespace std::complex_literals;

namespace Observable {

DensityContrast::DensityContrast(const Parameters& p, const Cosmology&)
    : sigma_x(p["Observables"]["DensityContrast"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Observables"]["DensityContrast"]["linear_convolution"]
                 .get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      t_prev(-1),
      ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
      delta(husimi ? 0 : N) {
    if (husimi) {
        auto& gaussian = ws.kernel_padded;
        std::iota(std::begin(gaussian), std::end(gaussian), -N_kernel / 2);
        gaussian =
            exp(-0.5 * dx * dx * gaussian * gaussian / (sigma_x * sigma_x));
        gaussian /= blaze::sum(gaussian);
    }
}

ObservableFunctor::ReturnType DensityContrast::compute(const SimState& state) {
    if (t_prev < state.tau) {
        t_prev = state.tau;
        if (husimi) {
            // Store delta inside convolution workspace
            ws.signal_padded = delta_from(state);
            discrete_convolution(ws);
            return ws.signal_padded;
        }
        // Store delta in dedicated array
        delta = delta_from(state);
    }

    return delta;
}

PhaseSpaceDistribution::PhaseSpaceDistribution(const Parameters& p,
                                               const Cosmology&)
    : sigma_x(p["Observables"]["PhaseSpaceDistribution"]["sigma_x"]),
      husimi(sigma_x > 0),
      linear(p["Observables"]["PhaseSpaceDistribution"]["linear_convolution"]),
      patch(p["Observables"]["PhaseSpaceDistribution"]["patch"]),
      N(p["Simulation"]["N"]),
      L(p["Simulation"]["L"]),
      N_kernel(2 * floor(5 * sigma_x / (L / N)) + 1),
      t_prev(-1),
      ws(0, 0, 0),
      wigner_f(0, 0),
      husimi_f(0, 0),
      idx(0),
      iaf(0, 0),
      c2c(nullptr) {
    if (husimi) {
        const double dx = L / N;
        const double dk = 2 * M_PI / L;

        // 0-based indices
        int idx_start = std::abs(-L / 2 - patch[0][0]) / dx;
        int idx_end = std::abs(-L / 2 - patch[0][1]) / dx;
        int idk_start = std::abs(-dk * N / 2 - patch[1][0]) / dk;
        int idk_end = std::abs(-dk * N / 2 - patch[1][1]) / dk;

        if (patch[0][0] >= patch[0][1] || patch[0][0] < -L / 2 ||
            L / 2 < patch[0][1]) {
            std::cout
                << WARNINGTAG(
                       "Wrong x interval in patch...Reset to entire domain")
                << std::endl;
            idx_start = 0;
            idx_end = N;
        }
        if (patch[1][0] >= patch[1][1] || patch[1][0] < -dk * N / 2 ||
            dk * (N - 1) / 2 < patch[1][1]) {
            std::cout
                << WARNINGTAG(
                       "Wrong k interval in patch...Reset to entire domain")
                << std::endl;
            idk_start = 0;
            idk_end = N;
        }

        const int Nx = idx_end - idx_start;
        const int Nk = idk_end - idk_start;

        idx.resize(Nx);
        std::iota(idx.begin(), idx.end(), -N / 2 + idx_start);
        idk.resize(Nk);
        std::iota(idk.begin(), idk.end(), -N / 2 + idk_start);
        husimi_f.resize(Nk, Nx);
        ws = convolution_ws<std::complex<double>>(linear, Nx, N_kernel);

        auto& gaussian = ws.kernel_padded;
        std::iota(std::begin(gaussian), std::end(gaussian), -N_kernel / 2);
        gaussian =
            exp(-dx * dx * gaussian * gaussian / (4 * sigma_x * sigma_x));
        gaussian /= blaze::sum(gaussian);
    } else {
        iaf.resize(N, N);
        wigner_f.resize(N, N);
        auto iaf_ptr = reinterpret_cast<fftw_complex*>(iaf.data());
        c2c.reset(fftw_plan_many_dft(
            1, &N, N, iaf_ptr, nullptr, 1, iaf.spacing(), iaf_ptr, nullptr, 1,
            iaf.spacing(), FFTW_FORWARD, FFTW_ESTIMATE));
    }
}

void PhaseSpaceDistribution::wigner_distribution(const SimState& state) {
    using namespace blaze;

    // A C++ version of Matlabs TFTB Toolkit Wigner-Ville-TFR

    auto negative_v_iaf = submatrix(iaf, N - N / 2, 0, N / 2, N);
    auto positive_v_iaf = submatrix(iaf, 0, 0, N - N / 2, N);

    auto negative_v_f = submatrix(wigner_f, 0, 0, N / 2, N);
    auto positive_v_f = submatrix(wigner_f, N - N / 2, 0, N - N / 2, N);

    auto& psi = state.psi;
    for (int icol = 0; icol < N; ++icol) {
        auto iaf_i = column(iaf, icol);
        // Caution: Round ties to next even integer value
        int ti = std::nearbyint(icol);
        int taumax = std::min(
            {ti, N - ti - 1, static_cast<int>(std::nearbyint(N / 2) - 1)});
        int total = 2 * taumax + 1;
        auto lag_plus = elements(
            psi, [&](int tau) { return icol + (tau - taumax); }, total);
        auto lag_minus = elements(
            psi, [&](int tau) { return icol - (tau - taumax); }, total);
        auto index = elements(
            iaf_i, [&](int tau) { return (N + tau - taumax) % N; }, total);
        index = lag_plus * conj(lag_minus);

        if (int tau = std::nearbyint(N / 2); ti <= N - tau && ti >= tau + 1)
            iaf(tau, icol) = 0.5 * (psi[ti + tau] * conj(psi[ti - tau]) +
                                    psi[ti - tau] * conj(psi[ti + tau]));
    }
    fftw_execute(c2c.get());

    negative_v_f = real(negative_v_iaf);
    positive_v_f = real(positive_v_iaf);
}

void PhaseSpaceDistribution::husimi_distribution(const SimState& state) {
    using namespace blaze;

    // Compute phase matrix only once
    static auto phase_T = exp(-1.0i * (2.0 * M_PI / N) * idx * trans(idk));

    // Construct mixed state distribution pure state by pure state
    auto& psi = state.psi;
    auto modulated_psi =
        phase_T %
        expand(subvector(psi, N / 2 + idx[0], idx.size()), idk.size());
    for (int k = 0; k < idk.size(); ++k) {
        ws.signal_padded = column(modulated_psi, k);
        discrete_convolution(ws);
        row(husimi_f, k) =
            trans(real(conj(ws.signal_padded) * ws.signal_padded));
    }
}

ObservableFunctor::ReturnType PhaseSpaceDistribution::compute(
    const SimState& state) {
    if (husimi) {
        if (t_prev < state.tau) {
            t_prev = state.tau;
            husimi_f = 0;
            husimi_distribution(state);
        }
        return husimi_f;
    }
    if (t_prev < state.tau) {
        t_prev = state.tau;
        wigner_f = 0;
        wigner_distribution(state);
    }
    return wigner_f;
}

Potential::Potential(const Parameters& p, const Cosmology&){};

inline ObservableFunctor::ReturnType Potential::compute(const SimState& state) {
    return state.V;
}

WaveFunction::WaveFunction(const Parameters& p, const Cosmology&){};

inline ObservableFunctor::ReturnType WaveFunction::compute(
    const SimState& state) {
    return state.psi;
}

Energy::Energy(const Parameters& p, const Cosmology& cosmo_)
    : cosmo(cosmo_),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      grad(N, N + 4),
      energies(4, 0),
      x(N) {
    grad.reserve(4 * N);

    // 5 point stencil for first derivative
    for (int i = 0; i < N; ++i) {
        grad.append(i, i, 1.0 / (12 * dx));
        grad.append(i, i + 1, -2.0 / (3 * dx));
        grad.append(i, i + 3, 2.0 / (3 * dx));
        grad.append(i, i + 4, -1.0 / (12 * dx));
        grad.finalize(i);
    }
    std::iota(x.begin(), x.end(), 0);
    x *= dx;
};

ObservableFunctor::ReturnType Energy::compute(const SimState& state) {
    using namespace blaze;
    auto cyclic_extension = [N = N](int i) { return (N - 2 + i) % N; };

    const double a = cosmo.a_of_tau(state.tau);

    // Kinetic Energy via trapezodial rule
    double& E_kinetic = energies[0];
    auto psi = elements(state.psi, cyclic_extension, N + 4);
    auto nabla_psi = grad * psi;
    E_kinetic = 0.5 * dx / (a * a) * sum(real(conj(nabla_psi) * nabla_psi));

    // Potential Energy via trapezodial rule
    double& E_pot = energies[1];
    E_pot = 0.5 * dx / a * sum(state.V * (real(conj(state.psi) * state.psi)));

    // Total Energy
    double& E_tot = energies[2];
    E_tot = E_kinetic + E_pot;

    // Assuming V to be arbitrary (i.e. potentially non homogeneous) we compute
    // the virial in its general form, i.e. <x grad V>
    double& x_grad_V = energies[3];
    auto V = elements(state.V, cyclic_extension, N + 4);
    auto nabla_V = grad * V;
    x_grad_V = dx / a * sum(x * nabla_V * (real(conj(state.psi) * state.psi)));

    return energies;
}

ParticleFlux::ParticleFlux(const Parameters& p, const Cosmology&)
    : N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      grad(N, N + 4),
      flux(N) {
    grad.reserve(4 * N);

    // 5 point stencil for first derivative
    for (int i = 0; i < N; ++i) {
        grad.append(i, i, 1.0 / (12 * dx));
        grad.append(i, i + 1, -2.0 / (3 * dx));
        grad.append(i, i + 3, 2.0 / (3 * dx));
        grad.append(i, i + 4, -1.0 / (12 * dx));
        grad.finalize(i);
    }
};

ObservableFunctor::ReturnType ParticleFlux::compute(const SimState& state) {
    using namespace blaze;
    auto cyclic_extension = [N = N](int i) { return (N - 2 + i) % N; };

    auto& psi = state.psi;
    auto psi_ext = elements(psi, cyclic_extension, N + 4);
    auto nabla_psi = grad * psi_ext;
    flux = real(0.5 * std::complex<double>{0, 1} *
                (psi * conj(nabla_psi) - conj(psi) * nabla_psi));

    return flux;
}

}  // namespace Observable

