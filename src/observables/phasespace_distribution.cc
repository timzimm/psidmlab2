#include <observables_common.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "logging.h"

namespace Observable {

class PhaseSpaceDistribution : public ObservableFunctor {
    const double sigma_x;  // spatial smoothing scale
    const bool husimi;     // compute husimi?
    const bool linear;     // do linear convolution?
    const Domain box;
    double t_prev;  // timestamp of last cached observable
    int idx_start;
    convolution_ws<std::complex<double>> ws;
    DynamicMatrix<double, columnMajor> wigner_f;  // cached wigner
    DynamicMatrix<double> husimi_f;               // cached husimi
    DynamicVector<double> x;
    DynamicVector<double> k;
    DynamicMatrix<std::complex<double>, columnMajor>
        iaf;            // instantaneous autocorrelation function
    fftw_plan_ptr c2c;  // IAF -> wigner transform (complex)

    // UNTESTED
    void wigner_distribution(const SimState& state) {
        // A C++ version of Matlabs TFTB Toolkit Wigner-Ville-TFR

        const int N = box.N;
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
        // Occasionally zeros crop up in f. Imagine a exponentially
        // decaying profile that drops below ~10^{-308} (smallest normalized
        // number). An exact 0, however, can lead to NaN effects in the derived
        // observables. We add the smallest possible offset to eliminate the
        // problem.
        const double eps = std::numeric_limits<double>::min();

        // Select part of psi that lives inside the x-patch and expand it in
        // k-direction
        auto psi = expand(subvector(state.psi, idx_start, x.size()), k.size());
        auto modulated_psi = psi % exp(-1.0i * x * trans(k));

        // columnwise gaussian smoothing + transpose into result matrix + offset
        for (int K = 0; K < k.size(); ++K) {
            ws.signal_padded = column(modulated_psi, K);
            discrete_convolution(ws);
            row(husimi_f, K) =
                eps + trans(real(conj(ws.signal_padded) * ws.signal_padded));
        }
    }

   public:
    PhaseSpaceDistribution(const Parameters& p, const Cosmology&)
        : sigma_x{p["Observables"]["PhaseSpaceDistribution"]["sigma_x"]},
          husimi(sigma_x > 0),
          linear{
              p["Observables"]["PhaseSpaceDistribution"]["linear_convolution"]},
          box(p),
          // Choose N_kernel such that it is 10 times the smoothing scale
          t_prev(-1),
          idx_start(0),
          ws(0, 0, 0),
          wigner_f(0, 0),
          husimi_f(0, 0),
          x(0),
          k(0),
          iaf(0, 0),
          c2c(nullptr) {
        if (husimi) {
            std::array<std::array<double, 2ul>, 2ul> patch;
            // Load patch information from config file.
            try {
                patch = p["Observables"]["PhaseSpaceDistribution"]["patch"];
            } catch (...) {
                std::cout << WARNINGTAG(
                                 "No valid patch specified. "
                                 "Use entire grid for Husimi function. "
                                 "This costs memory!")
                          << std::endl;
                // Default patch is entire grid
                patch = {{{box.xmin, box.xmax}, {box.kmin, box.kmax}}};
            }

            // Check if x-interval is...
            //(i)   in increasing order
            //(ii)  does not exceed the x-grid in negative direction
            //(iii) does not exceed the x-grid in positive direction
            //
            // In this context 'exceed' means the difference is larger than 100
            // times the machine epsilon
            const double threshold =
                100 * std::numeric_limits<double>::epsilon();
            if (patch[0][0] >= patch[0][1] ||
                box.xmin - patch[0][0] > threshold ||
                patch[0][1] - box.xmax > threshold) {
                std::cout
                    << WARNINGTAG(
                           "Wrong x interval in patch...Reset to entire domain")
                    << std::endl;
                patch[0][0] = box.xmin;
                patch[0][1] = box.xmax;
            }
            // Check if k-interval is...
            //(i)   in increasing order
            //(ii)  does not exceed the k-grid in negative direction
            //(iii) does not exceed the k-grid in positive direction
            //
            // In this context 'exceed' means the difference is larger than 100
            // times the machine epsilon

            if (patch[1][0] >= patch[1][1] ||
                box.kmin - patch[1][0] > threshold ||
                patch[1][1] - box.kmax > threshold) {
                std::cout
                    << WARNINGTAG(
                           "Wrong k interval in patch...Reset to entire domain")
                    << std::endl;
                patch[1][0] = box.kmin;
                patch[1][1] = box.kmax;
            }

            // index offset relative to the common lower boundary
            idx_start = static_cast<int>(
                std::abs(box.xmin - patch[0][0]) / box.dx + 0.5);
            int idk_start = static_cast<int>(
                std::abs(box.kmin - patch[1][0]) / box.dk + 0.5);
            int idx_end = static_cast<int>(
                std::abs(box.xmin - patch[0][1]) / box.dx + 0.5);
            int idk_end = static_cast<int>(
                std::abs(box.kmin - patch[1][1]) / box.dk + 0.5);

            const int Nx = idx_end - idx_start + 1;
            const int Nk = idk_end - idk_start + 1;

            x = box.xmin + box.dx * linspace(Nx, idx_start, idx_end);
            k = box.kmin + box.dk * linspace(Nk, idk_start, idk_end);

            husimi_f.resize(Nk, Nx);

            // Store gaussian kernel at most (floor) inside
            // [-8 sigma_x, +8 sigma_x] and make sure that N_kernel is odd
            // for symmetry of the kernel around x=0.
            int N_kernel = floor(16 * sigma_x / box.dx);
            N_kernel = N_kernel % 2 ? N_kernel : N_kernel + 1;
            if (N_kernel > box.N) {
                std::cerr << ERRORTAG("Smoothing scale is too large")
                          << std::endl;
                exit(1);
            }
            double L_kernel = (N_kernel - 1) * box.dx;
            auto xk = linspace(N_kernel, -L_kernel / 2, L_kernel / 2);

            ws = convolution_ws<std::complex<double>>(linear, Nx, N_kernel);

            auto& gaussian = ws.kernel_padded;
            // Gaussian smoothing kernel accoriding to fH defintion. We multiply
            // with dx to approximate the continuous convolution
            gaussian = box.dx / std::pow(2 * M_PI * sigma_x * sigma_x, 0.25) *
                       exp(-xk * xk / (4 * sigma_x * sigma_x));
        } else {
            iaf.resize(box.N, box.N);
            wigner_f.resize(box.N, box.N);
            auto iaf_ptr = reinterpret_cast<fftw_complex*>(iaf.data());
            c2c.reset(fftw_plan_many_dft(
                1, &box.N, box.N, iaf_ptr, nullptr, 1, iaf.spacing(), iaf_ptr,
                nullptr, 1, iaf.spacing(), FFTW_FORWARD, FFTW_ESTIMATE));
        }
    }

    ReturnType compute(
        const SimState& state,
        std::unordered_map<std::string, std::unique_ptr<ObservableFunctor>>&
            obs) {
        if (husimi) {
            if (t_prev < state.tau) {
                t_prev = state.tau;
                husimi_f = 0;
                husimi_distribution(state);
            }
            return husimi_f;
        } else {
            if (t_prev < state.tau) {
                t_prev = state.tau;
                wigner_f = 0;
                wigner_distribution(state);
            }
            return wigner_f;
        }
    }

    REGISTER(PhaseSpaceDistribution)
};

}  // namespace Observable
