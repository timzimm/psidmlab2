#include <observables_common.h>

namespace Observable {

class PhaseSpaceDistribution : public ObservableFunctor {
    using interval_t = std::array<double, 2ul>;
    using patch_t = std::array<interval_t, 2ul>;

    const double sigma_x;  // spatial smoothing scale
    const bool husimi;     // compute husimi?
    const bool linear;     // do linear convolution?
    const patch_t patch;
    const int N;         // number of spatial gridpoints
    const double L;      // spatial resolution
    const int N_kernel;  // symmetric 5-sigma_x interval in points
    double t_prev;       // timestamp of last cached observable
    convolution_ws<std::complex<double>> ws;
    DynamicMatrix<double, columnMajor> wigner_f;  // cached wigner
    DynamicMatrix<double> husimi_f;               // cached husimi
    DynamicVector<double> idx;
    DynamicVector<double> idk;
    DynamicMatrix<std::complex<double>, columnMajor>
        iaf;            // instantaneous autocorrelation function
    fftw_plan_ptr c2c;  // IAF -> wigner transform (complex)

    void wigner_distribution(const SimState& state) {
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

    void husimi_distribution(const SimState& state) {
        // Compute phase matrix only once
        static auto phase_T = exp(-1.0i * (2.0 * M_PI / N) * idx * trans(idk));

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

   public:
    PhaseSpaceDistribution(const Parameters& p, const Cosmology&)
        : sigma_x{p["Observables"]["PhaseSpaceDistribution"]["sigma_x"]},
          husimi(sigma_x > 0),
          linear{
              p["Observables"]["PhaseSpaceDistribution"]["linear_convolution"]},
          patch{p["Observables"]["PhaseSpaceDistribution"]["patch"]},
          N{p["Simulation"]["N"]},
          L{p["Simulation"]["L"]},
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
            gaussian /= sum(gaussian);
        } else {
            iaf.resize(N, N);
            wigner_f.resize(N, N);
            auto iaf_ptr = reinterpret_cast<fftw_complex*>(iaf.data());
            c2c.reset(fftw_plan_many_dft(
                1, &N, N, iaf_ptr, nullptr, 1, iaf.spacing(), iaf_ptr, nullptr,
                1, iaf.spacing(), FFTW_FORWARD, FFTW_ESTIMATE));
        }
    }

    ReturnType compute(
        const SimState& state,
        const std::unordered_map<std::string,
                                 std::unique_ptr<ObservableFunctor>>& obs) {
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

    REGISTER(PhaseSpaceDistribution)
};

}  // namespace Observable
