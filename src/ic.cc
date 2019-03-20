#include "ic.h"
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include "fftw_alloc.h"

ICGenerator::ICGenerator(const Parameters& param)
    : type(param.ic),
      N(param.N),
      dx(param.dx),
      L(param.L),
      M(param.M),
      rel_threshold(param.ev_thr),
      source_name(param.ic_source_file),
      ic_file(source_name) {
    // Altpugh the input file format is generic for rho and powerspectrum
    // type initial conditions, we leave the exact way of parsing the data to
    // the generator routines to allow for specialized malloc's such as done by
    // FFTW3. Here, we only discard the file header and count lines
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');

    // Number of lines
    data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                        std::istreambuf_iterator<char>(), '\n');
    data_N = (data_N + 1);
    ic_file.seekg(0);
    ic_file.ignore(std::numeric_limits<int>::max(), '\n');
}

// Dispatches the correct generator routine depending on
// Parameters::IC parameter at runtime
void ICGenerator::generate(SimState& state) const {
    if (type == Parameters::IC::Density)
        psi_from_rho(state);
    else
        psi_from_power(state);
}

void ICGenerator::psi_from_rho(SimState& state) const {
    // Read in density data from file
    std::vector<double, FFTWAlloc<double>> rho(data_N);
    std::vector<std::complex<double>, FFTWAlloc<std::complex<double>>> rho_fft(
        data_N / 2 + 1);
    for (int i = 0; i < data_N; ++i) {
        ic_file >> rho[i] >> rho[i];
        if (!ic_file) break;
    }

    // Compute FFT to obtain coefficients for the expansion of rho in
    // harmonic base functions
    auto plan = fftw_plan_dft_r2c_1d(
        data_N, rho.data(), reinterpret_cast<fftw_complex*>(rho_fft.data()),
        FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // The modulus of the Fourier coefficients are the eigenvalues of the
    // integral operator eigenvalue problem to solve. We only care about
    // dominant modes in the statistic mixture, hence we order all coefficients
    // by magnitude. In fact, we also need information about the mode number n
    // to sample from the correct harmonics later on. Therefore, we compute the
    // permution index array that yields a sorted coefficient list.

    std::vector<int> mode(rho_fft.size());
    std::iota(mode.begin(), mode.end(), 0);

    // Sort in descending order of coefficient modulus
    std::sort(mode.begin(), mode.end(), [&](const int a, const int b) {
        return (norm(rho_fft[a]) < norm(rho_fft[b]));
    });

    int M_true;
    // Deduce operation type...
    // (i) rel_threshold <= 0 -> min(M, #fourier_modes)
    if (rel_threshold <= 0) {
        std::cout << INFOTAG("Determine No. of wavefunctions by fiven M")
                  << std::endl;
        M_true = std::min(M, static_cast<int>(rho_fft.size()));
    }

    // (ii) rel_threshold > 0 -> min(M s.t. all ev >
    // ev_threshold,#fourier_modes)
    else {
        std::cout << INFOTAG(
                         "Determine No. of wavefunctions by given threshold")
                  << std::endl;
        // min(..) is already taken care of.
        M_true = std::count_if(mode.begin(), mode.end(), [&](const int m) {
            return rel_threshold * norm(rho_fft[0]) < norm(rho_fft[m]);
        });
    }

    std::cout << INFOTAG("Construct" << M_true << " wavefunctions")
              << std::endl;

    // Construct initial wave functions
    // For each mode there are two equally dominant solutions namely
    // lambda(+) and lambda(i)
    for (int m = 0; m < M_true; m += 2) {
        int n = mode[m / 2];
        std::complex<double> c = rho_fft[n];
        double a = 2 * c.real();
        double b = -2 * c.imag();

        double lambda = abs(c);
        double normal = 1 / std::sqrt((a + lambda) * (a + lambda) + b * b);
        // Append (+)-solution
        for (int i = 0; i < N; ++i) {
            state.psis.push_back(
                normal * ((a + lambda) * std::cos(2 * M_PI * n * dx * i / L) +
                          b * std::sin(2 * M_PI * n * dx * i / L)));
        }

        lambda *= -1.0;
        normal = 1 / std::sqrt((a + lambda) * (a + lambda) + b * b);
        // Append (-)-solution
        for (int i = 0; i < N; ++i) {
            state.psis.push_back(
                normal * ((a + lambda) * std::cos(2 * M_PI * n * dx * i / L) +
                          b * std::sin(2 * M_PI * n * dx * i / L)));
        }
    }
}
void ICGenerator::psi_from_power(SimState& state) const {}
