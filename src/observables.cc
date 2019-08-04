#include "observables.h"
#include "blaze_utils.h"
#include "parameters.h"
#include "state.h"

#include <blaze/math/Columns.h>
#include <blaze/math/Elements.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Submatrix.h>
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

PhaseSpaceDistribution::PhaseSpaceDistribution(const Parameters& p)
    : sigma_x(std::sqrt(2) * p["Observables"]["sigma_x"].get<double>()),
      husimi(sigma_x > 0),
      linear(p["Observables"]["linear_convolution"].get<bool>()),
      N(p["Simulation"]["N"].get<int>()),
      dx(p["Simulation"]["L"].get<double>() / N),
      N_kernel(2 * floor(5 * sigma_x / dx) + 1),
      t_prev(-1),
      ws(husimi ? linear : 0, husimi ? N : 0, husimi ? N_kernel : 0),
      wigner_f(husimi ? 0 : N, husimi ? 0 : N),
      husimi_f(husimi ? N : 0, husimi ? N : 0),
      idx(husimi ? N : 0),
      iaf(husimi ? 0 : N, husimi ? 0 : N),
      c2c(nullptr) {
    if (husimi) {
        std::iota(idx.begin(), idx.end(), 0);
        auto& gaussian = ws.kernel_padded;
        std::iota(std::begin(gaussian), std::end(gaussian), -N_kernel / 2);
        gaussian =
            exp(-0.5 * dx * dx * gaussian * gaussian / (sigma_x * sigma_x));
        gaussian /= blaze::sum(gaussian);
    } else {
        auto iaf_ptr = reinterpret_cast<fftw_complex*>(iaf.data());
        c2c = fftw_plan_many_dft(1, &N, N, iaf_ptr, nullptr, 1, iaf.spacing(),
                                 iaf_ptr, nullptr, 1, iaf.spacing(),
                                 FFTW_FORWARD, FFTW_ESTIMATE);
    }
}

void PhaseSpaceDistribution::wigner_distribution(const SimState& state) {
    using namespace blaze;

    // A C++ version of Matlabs TFTB Toolkit Wigner-Ville-TFR

    auto negative_v_iaf = submatrix(iaf, N - N / 2, 0, N / 2, N);
    auto positive_v_iaf = submatrix(iaf, 0, 0, N - N / 2, N);

    auto negative_v_f = submatrix(wigner_f, 0, 0, N / 2, N);
    auto positive_v_f = submatrix(wigner_f, N - N / 2, 0, N - N / 2, N);

    for (int m = 0; m < state.M; ++m) {
        auto psi = column(state.psis, m);
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
        fftw_execute(c2c);

        negative_v_f += state.lambda[m] * real(negative_v_iaf);
        positive_v_f += state.lambda[m] * real(positive_v_iaf);
    }
}

void PhaseSpaceDistribution::husimi_distribution(const SimState& state) {
    using namespace blaze;

    // Compute phase matrix only once
    static auto phase_T = exp(idx * (-M_PI + 2 * M_PI / N * trans(idx)) *
                              std::complex<double>(0, -1));

    // Construct mixed state distribution pure state by pure state
    for (int m = 0; m < state.M; ++m) {
        auto psi = column(state.psis, m);
        auto modulated_psi = phase_T % expand(psi, N);
        for (int k = 0; k < N; ++k) {
            ws.signal_padded = column(modulated_psi, k);
            discrete_convolution(ws);
            row(husimi_f, k) +=
                state.lambda[m] *
                trans(real(conj(ws.signal_padded) * ws.signal_padded));
        }
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

PhaseSpaceDistribution::~PhaseSpaceDistribution() { fftw_destroy_plan(c2c); }

Potential::Potential(const Parameters& p)
    : pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      potential(p["Simulation"]["N"].get<int>()){};

inline ObservableFunctor::ReturnType Potential::compute(const SimState& state) {
    // Depending on the integrator, state.V might be in an intermediate step and
    // only state.psis is correct. Thus, we recalculate the potential.
    pot->solve(potential, delta_from(state));
    return potential;
}

WaveFunction::WaveFunction(const Parameters& p){};

inline ObservableFunctor::ReturnType WaveFunction::compute(
    const SimState& state) {
    return state.psis;
}

Energy::Energy(const Parameters& p)
    : N(p["Simulation"]["N"].get<int>()),
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

    // Kinetic Energy via trapezodial rule
    double& E_kinetic = energies[0];
    auto psi = rows(state.psis, cyclic_extension, N + 4);
    auto nabla_psi = grad * psi;
    E_kinetic = 0.5 * dx / (state.a * state.a) *
                sum(real(conj(nabla_psi) % nabla_psi) * state.lambda);

    // Potential Energy via trapezodial rule
    double& E_pot = energies[1];
    E_pot = 0.5 * dx / state.a *
            sum(state.V * (real(conj(state.psis) % state.psis) * state.lambda));

    // Total Energy
    double& E_tot = energies[2];
    E_tot = E_kinetic + E_pot;

    // Assuming V to be arbitrary (i.e. potentially non homogeneous) we compute
    // the virial in its general form, i.e. <x grad V>
    double& x_grad_V = energies[3];
    auto V = elements(state.V, cyclic_extension, N + 4);
    auto nabla_V = grad * V;
    x_grad_V =
        dx / state.a *
        sum(x * nabla_V * (real(conj(state.psis) % state.psis) * state.lambda));

    return energies;
}

ParticleFlux::ParticleFlux(const Parameters& p)
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

// Currently only M=1
ObservableFunctor::ReturnType ParticleFlux::compute(const SimState& state) {
    using namespace blaze;
    auto cyclic_extension = [N = N](int i) { return (N - 2 + i) % N; };

    auto psi = column(state.psis, 0);
    auto psi_ext = elements(psi, cyclic_extension, N + 4);
    auto nabla_psi = grad * psi_ext;
    flux = real(0.5 * std::complex<double>{0, 1} *
                (psi * conj(nabla_psi) - conj(psi) * nabla_psi));

    return flux;
}

}  // namespace Observable

