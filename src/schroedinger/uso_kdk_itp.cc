#include "schroedinger/uso_kdk_itp.h"
#include "cosmology.h"
#include "io.h"
#include "parameters.h"

#include <algorithm>
#include <fstream>

using namespace blaze;
using namespace std::complex_literals;

namespace Schroedinger {

USO_KDK_ITP::USO_KDK_ITP(const Parameters& p, const SimState& state,
                         const Cosmology& cosmo_)
    : DefaultDriver{p},
      cosmo{cosmo_},
      pot{PotentialMethod::make(p["Simulation"]["potential"].get<std::string>(),
                                p)},
      pot_external(0),
      N{p["Simulation"]["N"].get<int>()},
      L{p["Simulation"]["L"].get<double>()},
      norm0{sqrt(L / N *
                 sum(real(conj(state.psis) % state.psis) * state.lambda))},
      k_squared(N),
      kick(N, state.M),
      dt_last{-1} {
    // Compute initial normalisation
    // FFTW reorders frequencies. The upper half starts at the most negative
    // frequency and increases for increasing index.
    //          f0 f1 ... fN/2-1, f-N/2 ... f-1
    std::iota(k_squared.begin(), k_squared.end(), -N / 2);
    std::rotate(k_squared.begin(), k_squared.begin() + N / 2, k_squared.end());
    k_squared = 4 * M_PI * M_PI / (L * L) * k_squared * k_squared;

    // Check if an external potential is required
    std::ifstream pot_file{
        p["Simulation"]["stepper"]["external_potential"].get<std::string>()};
    if (!pot_file) return;
    std::cout << INFOTAG("External potential loaded") << std::endl;

    // Determine # rows in file
    double dataN = std::count(std::istreambuf_iterator<char>(pot_file),
                              std::istreambuf_iterator<char>(), '\n');
    if (dataN != N) {
        std::cout << ERRORTAG("N from external potential differs") << std::endl;
        exit(1);
    }
    pot_external.resize(N);
    pot_file.clear();
    pot_file.seekg(0, std::ios::beg);
    fill_from_file(pot_file, pot_external);
}

void USO_KDK_ITP::step(SimState& state, const double dt) {
    CCM& psis = state.psis;
    const double a_half = cosmo.a_of_tau(state.tau + dt / 2);

    // Set kick operator once for the entire time step
    if (dt - dt_last != 0) {
        kick = expand(exp(-0.5 * 0.5 * k_squared * dt), state.M);
    }

    state.transform(SimState::Representation::Momentum);
    psis = kick % psis;

    state.transform(SimState::Representation::Position);
    // Recompute potential
    pot->solve(state);

    // Set drift operator with intermediate potential
    if (pot_external.size() == N) state.V += pot_external;
    auto drift = expand(exp(-1.0 * a_half * state.V * dt), state.M);
    psis = drift % psis;

    state.transform(SimState::Representation::Momentum);
    psis = kick % psis;
    state.transform(SimState::Representation::Position);

    double norm2 =
        L / N * sum(real(conj(state.psis) % state.psis) * state.lambda);

    // Renormalize to box size
    column(state.psis, 0) /= sqrt(norm2) / norm0;

    state.tau += dt;
    dt_last = dt;
    state.n += 1;
}

double USO_KDK_ITP::next_dt(const SimState& state) const {
    return std::min(
        {L * L / (N * N * M_PI), M_PI / (cosmo.a_of_tau(state.tau) *
                                         max(abs(state.V + pot_external)))});
}

}  // namespace Schroedinger
