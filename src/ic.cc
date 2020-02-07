#include "ic.h"
#include "cosmology.h"
#include "fftw.h"
#include "gsl_integration.h"
#include "interfaces.h"
#include "io.h"
#include "parameters.h"
#include "state.h"

#include <algorithm>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <fstream>
#include <numeric>
#include <random>
#include <vector>
// Special Functions for integration
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_trig.h>

using namespace blaze;
using namespace boost::math;

void set_phase(SimState&, const Cosmology&, const Parameters&);
void psi_from_file(DynamicVector<std::complex<double>>&, std::ifstream&);
void psi_polar_from_file(DynamicVector<std::complex<double>>&, std::ifstream&);
void delta_from_power(SimState&, const Cosmology&, std::ifstream& k_Pk);

template <typename Functor>
auto reduced_powerspectrum(const int dim, const Functor& P3d) {
    // Minimum k decade from which P is computed
    const double decmin = -4;
    // Maximum k decade from which P is computed
    const double decmax = 2;
    // Samples per decade
    const double Nperlog = 100;
    // Total number of samples
    const double N = Nperlog * (decmax - decmin);
    // Sampling increment
    const double dlogk = (decmax - decmin) / (N - 1);

    double logk = decmin;
    double k_ortho = std::pow(10, decmin);

    // Maximum number of subintervals used in integration
    const int n = 200;
    // Starting relative accuracy + increment
    const double eps = 1e-7;
    const double eps_increment = 1e-7;

    // Sampled values of correlation function.
    blaze::DynamicVector<double> Ps(N, 0);

    auto ws = gsl_integration_workspace_alloc(n);

    auto integration_loop = [&](auto f) {
        for (auto& P : Ps) {
            double eps_local = eps;
            double err;
            while (gsl_integration_qagiu(f, k_ortho, 0, eps_local, n, ws, &P,
                                         &err) == GSL_EROUND) {
                eps_local += eps_increment;
            }
            logk += dlogk;
            k_ortho = std::pow(10, logk);
        }
    };

    if (dim == 1) {
        auto f = make_gsl_function([&](double k) { return k * P3d(k); });
        integration_loop(f);
        Ps *= 1.0 / (2 * M_PI);
    } else if (dim == 2) {
        auto f = make_gsl_function([&](double k) {
            return k * P3d(k) / sqrt(k * k + k_ortho * k_ortho);
        });
        integration_loop(f);
        Ps *= 1.0 / M_PI;
    }

    gsl_integration_workspace_free(ws);
    auto P = [dim = dim,
              P = boost::math::cubic_b_spline<double>(Ps.begin(), Ps.end(),
                                                      decmin, dlogk),
              kmin = std::pow(10, decmin),
              kmax = std::pow(10, decmax)](double k) {
        // Power law extrapolation before kmin. Assum P(k) = A*k^0
        if (k < kmin) {
            return P(std::log10(kmin));
        }
        // Power law extrapolation past kmax. Assume P(k) = A*k^-dim
        if (k > kmax) {
            return P(std::log10(kmax)) * std::pow(kmax, dim) *
                   std::pow(k, -dim);
        }
        return P(std::log10(k));
    };
    return P;
}

template <typename Functor>
auto correlation_function(const int dim, const double R, const Functor& P) {
    // Frequency after which the power spectrum takes on a power law
    const double kmax = 100;

    // Maximum number of subintervals used in integration
    const int n = 300;
    // Starting relative accuracy + increment
    const double eps = 1e-7;
    const double eps_increment = 1e-7;

    const double rmin = 1e-4;
    const double logrmin = std::floor(std::log10(rmin));

    const double rcut = 10;
    const double logrcut = std::log10(rcut);

    // Total number of samples
    const int Nlog = 100;
    const int N = 100;

    // Sampling increments
    const double dlogr = (logrcut - logrmin) / (Nlog - 1);
    const double dr = (R / 2 - rcut) / (N - 1);

    double logr = logrmin;
    double r = std::pow(10, logr);

    // Linearly sampled values of correlation function.
    blaze::DynamicVector<double> xis(N, 0);
    blaze::DynamicVector<double> xis_logr(N, 0);

    auto ws = gsl_integration_workspace_alloc(n);
    auto tab = gsl_integration_qawo_table_alloc(r, 0, GSL_INTEG_SINE, n);

    switch (dim) {
        case 1: {
            auto f = make_gsl_function([&](double k) { return P(k); });
            for (auto& xi : xis_logr) {
                double err;
                double eps_local = eps;
                gsl_integration_qawo_table_set(tab, r, kmax, GSL_INTEG_COSINE);
                while (gsl_integration_qawo(f, 0.0, 0, eps_local, n, ws, tab,
                                            &xi, &err) == GSL_EROUND) {
                    eps_local += eps_increment;
                }
                xi -= kmax * P(kmax) * gsl_sf_Ci(kmax * r);
                logr += dlogr;
                r = std::pow(10, logr);
            }
            r = rcut;
            for (auto& xi : xis) {
                double err;
                double eps_local = eps;
                gsl_integration_qawo_table_set(tab, r, kmax, GSL_INTEG_COSINE);
                while (gsl_integration_qawo(f, 0.0, 0, eps_local, n, ws, tab,
                                            &xi, &err) == GSL_EROUND) {
                    eps_local += eps_increment;
                }
                xi -= kmax * P(kmax) * gsl_sf_Ci(kmax * r);
                r += dr;
            }
            xis_logr *= 1.0 / M_PI;
            xis *= 1.0 / M_PI;
            break;
        }
        case 2: {
            auto f = make_gsl_function(
                [&](double k) { return k * gsl_sf_bessel_J0(k * r) * P(k); });

            const double xi_threshold = 1e-6;
            const double k_threshold = 200;
            for (auto& xi : xis_logr) {
                int i = 3;
                double err;
                double xi_part = 0.0;
                double lower = 0.0;
                double upper = gsl_sf_bessel_zero_J0(i);
                while (lower < k_threshold ||
                       abs(xi_part) < xi_threshold * abs(xi)) {
                    double eps_local = eps;
                    while (gsl_integration_qag(f, lower, upper, 0, eps_local, n,
                                               GSL_INTEG_GAUSS21, ws, &xi_part,
                                               &err) == GSL_EROUND) {
                        eps_local += eps_increment;
                    }
                    xi += xi_part;
                    lower = upper;
                    i += 3;
                    upper = gsl_sf_bessel_zero_J0(i);
                }
                logr += dlogr;
                r = std::pow(10, logr);
            }
            r = rcut;
            for (auto& xi : xis) {
                int i = 3;
                double err;
                double xi_part = 0.0;
                double lower = 0.0;
                double upper = gsl_sf_bessel_zero_J0(i);
                while (lower < k_threshold ||
                       abs(xi_part) < xi_threshold * abs(xi)) {
                    double eps_local = eps;
                    while (gsl_integration_qag(f, lower, upper, 0, eps_local, n,
                                               GSL_INTEG_GAUSS21, ws, &xi_part,
                                               &err) == GSL_EROUND) {
                        eps_local += eps_increment;
                    }
                    xi += xi_part;
                    lower = upper;
                    i += 3;
                    upper = gsl_sf_bessel_zero_J0(i);
                }
                r += dr;
            }
            xis_logr *= 1.0 / (2 * M_PI);
            xis *= 1.0 / (2 * M_PI);
            break;
        }
        case 3: {
            auto f = make_gsl_function([&](double k) { return P(k) * k / r; });

            for (auto& xi : xis_logr) {
                double err;
                double eps_local = eps;
                gsl_integration_qawo_table_set(tab, r, kmax, GSL_INTEG_SINE);
                while (gsl_integration_qawo(f, 0.0, 0, eps_local, n, ws, tab,
                                            &xi, &err) == GSL_EROUND) {
                    eps_local += eps_increment;
                }
                xi += std::pow(kmax, 3) * P(kmax) *
                      (gsl_sf_sinc(kmax * r / M_PI) - gsl_sf_Ci(kmax * r));
                logr += dlogr;
                r = std::pow(10, logr);
            }
            r = rcut;
            for (auto& xi : xis) {
                double err;
                double eps_local = eps;
                gsl_integration_qawo_table_set(tab, r, kmax, GSL_INTEG_SINE);
                while (gsl_integration_qawo(f, 0.0, 0, eps_local, n, ws, tab,
                                            &xi, &err) == GSL_EROUND) {
                    eps_local += eps_increment;
                }
                xi += std::pow(kmax, 3) * P(kmax) *
                      (gsl_sf_sinc(kmax * r / M_PI) - gsl_sf_Ci(kmax * r));
                r += dr;
            }
            xis_logr *= 1 / (2 * M_PI * M_PI);
            xis *= 1 / (2 * M_PI * M_PI);
            break;
        }
    }

    gsl_integration_workspace_free(ws);
    gsl_integration_qawo_table_free(tab);

    // Interpolator for correlation function
    auto xi = [Plin = cubic_b_spline<double>(xis.begin(), xis.end(), rcut, dr),
               Plog = cubic_b_spline<double>(xis_logr.begin(), xis_logr.end(),
                                             logrmin, dlogr),
               rcut = rcut](double r) {
        // Technically wrong, but doesn't matter here as xi is multiplied by
        // r in 2,3D later on in the integration.
        if (r == 0.0) {
            return 0.0;
        }
        if (r < rcut) {
            return Plog(std::log10(r));
        }
        return Plin(r);
    };
    return xi;
}

template <typename Functor>
auto convolved_powerspectrum(const int dim, const double R, const Functor& xi) {
    // Minimum k from which P is computed
    const double decmin = -4;
    // Maximum k from which P is computed
    const double decmax = 2;
    // Samples per decade
    const double Nperlog = 100;
    // Total number of samples
    const double N = Nperlog * (decmax - decmin);
    // Sampling increment
    const double dlogk = (decmax - decmin) / (N - 1);

    double logk = decmin;
    double k = std::pow(10, decmin);

    // Maximum number of subintervals used in integration
    const int n = 200;
    // Starting relative accuracy + increment
    const double eps = 1e-7;
    const double eps_increment = 1e-7;

    // Log sampled values of correlation function.
    blaze::DynamicVector<double> Ps(N, 0);

    auto ws = gsl_integration_workspace_alloc(n);
    auto tab = gsl_integration_qawo_table_alloc(k, R / 2, GSL_INTEG_SINE, n);

    auto fourier_transformation = [&](auto f, auto type) {
        for (auto& P : Ps) {
            double err;
            double eps_local = eps;
            gsl_integration_qawo_table_set(tab, k, R / 2, type);
            while (gsl_integration_qawo(f, 0.0, 0, eps_local, n, ws, tab, &P,
                                        &err) == GSL_EROUND) {
                eps_local += eps_increment;
            }
            logk += dlogk;
            k = std::pow(10, logk);
        }
    };
    // Compute box adjusted power spectrum depending on space dimension dim.
    switch (dim) {
        case 1: {
            auto f = make_gsl_function([&](double r) { return xi(r); });
            fourier_transformation(f, GSL_INTEG_COSINE);
            Ps *= 2;
            break;
        }
        case 2: {
            auto f = make_gsl_function(
                [&](double r) { return gsl_sf_bessel_J0(k * r) * xi(r) * r; });

            for (auto& P : Ps) {
                double err;
                double P_part = 0.0;
                double lower = 0.0;
                double upper = gsl_sf_bessel_zero_J0(3);
                for (size_t i = 2; lower < R / 2; i += 3) {
                    double eps_local = eps;
                    while (gsl_integration_qag(f, lower, upper, 0, eps_local, n,
                                               GSL_INTEG_GAUSS21, ws, &P_part,
                                               &err) == GSL_EROUND) {
                        eps_local += eps_increment;
                    }
                    P += P_part;
                    lower = upper;
                    upper = std::min(gsl_sf_bessel_zero_J0(i), R / 2);
                }
                logk += dlogk;
                k = std::pow(10, logk);
            }
            Ps *= 2 * M_PI;
            break;
        }
        case 3: {
            auto f = make_gsl_function([&](double r) { return xi(r) * r / k; });
            fourier_transformation(f, GSL_INTEG_SINE);
            Ps *= 4 * M_PI;
            break;
        }
    }
    gsl_integration_workspace_free(ws);
    gsl_integration_qawo_table_free(tab);

    auto P = [dim = dim,
              P = boost::math::cubic_b_spline<double>(Ps.begin(), Ps.end(),
                                                      decmin, dlogk),
              kmin = std::pow(10, decmin),
              kmax = std::pow(10, decmax)](double k) {
        // Power law extrapolation before kmin. Assum P(k) = A*k^0
        if (k < kmin) {
            return P(std::log10(kmin));
        }
        // Power law extrapolation past kmax. Assume P(k) = A*k^-dim
        if (k > kmax) {
            return P(std::log10(kmax)) * std::pow(kmax, dim) *
                   std::pow(k, -dim);
        }
        return P(std::log10(k));
    };

    return P;
}

void generate_ic(SimState& state, const Cosmology& cosmo, const Parameters& p) {
    auto type = static_cast<ICType>(p["Initial Conditions"]["ic_type"]);
    std::ifstream ic_file{
        p["Initial Conditions"]["source_file"].get<std::string>()};

    switch (type) {
        case ICType::RealImagPsi... ICType::PolarPsi: {
            size_t data_N = std::count(std::istreambuf_iterator<char>(ic_file),
                                       std::istreambuf_iterator<char>(), '\n');
            if (data_N != state.box.N_total) {
                std::cout << ERRORTAG("#lines in source_file(" << data_N
                                                               << ") != N")
                          << std::endl;
                exit(1);
            }
            ic_file.seekg(0);
        }
        default:
            break;
    };

    switch (type) {
        case ICType::RealImagPsi:
            std::cout << INFOTAG("Load Re(psi0)/Im(psi0) from file")
                      << std::flush;
            psi_from_file(state.psi, ic_file);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::PolarPsi:
            std::cout << INFOTAG("Load Psi0 from file") << std::flush;
            psi_polar_from_file(state.psi, ic_file);
            std::cout << " ... done" << std::endl;
            break;
        case ICType::Powerspectrum:
            std::cout << INFOTAG("Generate delta0 from P(k)") << std::flush;
            // TODO: Don't return delta in state.V
            delta_from_power(state, cosmo, ic_file);
            std::cout << " ... done" << std::endl;
            break;
    }

    // If psi comes from a file its velocity is already set. For powerspectra
    // initializing the velocity is a additional (potentially unwanted) step
    // Velocity Initialization
    if (type == ICType::Powerspectrum) {
        // state.V holds delta at this point.
        // Set psis modulus before we override delta. For dust-like IC
        // all that is left to do.
        state.psi = sqrt(1.0 + state.V);

        if (p["Initial Conditions"]["compute_velocity"].get<bool>()) {
            std::cout << INFOTAG("Initial Velocity Field from Poisson ")
                      << std::endl;
            set_phase(state, cosmo, p);
        }
    }
}

void set_phase(SimState& state, const Cosmology& cosmo, const Parameters& p) {
    auto& delta = state.V;
    auto& phase = state.V;
    auto& psi = state.psi;
    auto poisson = Interaction::make("Poisson::FFT", p);
    const double a_init = cosmo.a_of_tau(0);
    const double prefactor =
        -std::sqrt(2.0 / 3 * a_init / cosmo.omega_m(a_init));

    poisson->solve(phase, prefactor * delta);

    // Madelung Representation.
    psi *= exp(std::complex<double>(0, 1) * phase);
}

void psi_from_file(DynamicVector<std::complex<double>>& psi, std::ifstream& f) {
    fill_from_file(f, psi);
}

void psi_polar_from_file(DynamicVector<std::complex<double>>& psi,
                         std::ifstream& f) {
    fill_from_file(f, psi);
    psi = map(psi, [](std::complex<double> c) {
        return std::polar(c.real(), c.imag());
    });
}

// TODO: k/h, kh, k, what the hell?
void delta_from_power(SimState& state, const Cosmology& cosmo,
                      std::ifstream& k_Pk) {
    const size_t N = std::count(std::istreambuf_iterator<char>(k_Pk),
                                std::istreambuf_iterator<char>(), '\n');
    k_Pk.clear();
    k_Pk.seekg(0, std::ios::beg);

    DynamicVector<double> ks(N);
    DynamicVector<double> Pk(N);
    fill_from_file(k_Pk, ks, Pk);

    double logk0 = std::log10(ks[0]);
    double delta_logk = std::log10(ks[1]) - logk0;

    auto P_3d = [P = cubic_b_spline<double>(Pk.begin(), Pk.end(), logk0,
                                            delta_logk),
                 kmax = ks[N - 1], kmin = ks[0]](double k) {
        // Power law extrapolation before kmin. Assum P(k) = A*k
        if (k < kmin) {
            return P(std::log10(kmin)) / kmin * k;
        }
        // Power law extrapolation past kmax. Assume P(k) = A*k^-3
        if (k > kmax) {
            return P(std::log10(kmax)) * std::pow(kmax, 3) * std::pow(k, -3);
        }
        return P(std::log10(k));
    };

    const Domain box = state.box;
    const int Nx = box.Ns[0];
    const int Ny = box.Ns[1];
    const int Nz = box.Ns[2];

    auto P_dim = reduced_powerspectrum(box.dims, P_3d);

    const auto& Ls_pc = box.box_lengths_pc;
    const double L = sqrt(
        std::inner_product(Ls_pc.begin(), Ls_pc.end(), Ls_pc.begin(), 0.0));
    auto xi = correlation_function(box.dims, L / 1e6, P_dim);

    auto P = convolved_powerspectrum(box.dims, L, xi);

    const size_t N_modes = (box.N_total / Nx) * (Nx / 2 + 1);
    auto deltas = subvector(state.psi, 0, N_modes);
    auto kx = subvector(state.V, 0, Nx / 2 + 1);
    auto ky = subvector(state.V, Nx / 2 + 1, Ny);
    auto kz = subvector(state.V, Ny + Nx / 2 + 1, Nz);

    // By hermitian symmetry of the real to complex transform only half positive
    // frequencies are relevant, i.e.
    //                      k_0, k_1, ... ,k_N/2
    std::iota(kx.begin(), kx.end(), 0.0);

    // For all other dimensions we have to consider the entire frequency domain
    // since after transforming in x direction the coefficients are complex
    // again, hence no symmetry can be assumed.
    //         k_0, k_1, ... ,k_N/2, k_(-N/2+1), ... , k_-1
    // For dummy dimensions (Ni=1) we round up to 2 in order to get k=0 as only
    // mode
    auto next_even = [](int i) {
        return static_cast<int>(std::round(i / 2.0) * 2.0);
    };
    std::iota(ky.begin(), ky.end(), -1 * next_even(Ny) / 2 + 1);
    std::rotate(ky.begin(), ky.begin() + next_even(Ny) / 2 - 1, ky.end());
    std::iota(kz.begin(), kz.end(), -1 * next_even(Nz) / 2 + 1);
    std::rotate(kz.begin(), kz.begin() + next_even(Nz) / 2 - 1, kz.end());

    // Finally, since P(k) expects k to be in h * Mpc^(-1), we multiply by the
    // correct increment
    kx *= 2 * M_PI / Ls_pc[0] * 1e6;
    ky *= 2 * M_PI / Ls_pc[1] * 1e6;
    kz *= 2 * M_PI / Ls_pc[2] * 1e6;

    // The entire k-grid can than be represented as sum of Kronecker
    // products.
    auto ex = uniform(Nx / 2 + 1, 1.0);
    auto ey = uniform(Ny, 1.0);
    auto ez = uniform(Ny, 1.0);

    auto k_mag =
        sqrt(kron(kron(ez, ey), kx * kx) + kron(kron(ez, ky * ky), ex) +
             kron(kron(kz * kz, ey), ex));

    // Independent number generators for modulus and phase.
    std::random_device rd;
    auto u1 =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(rd()));
    auto u2 =
        std::bind(std::uniform_real_distribution<>{0, 1}, std::mt19937(rd()));

    // Rescale to initial time according to linear theory encapsulated in D(a)
    const double a_init = cosmo.a_of_tau(0);
    const double D_init_0 = cosmo.D(a_init) / cosmo.D(1);

    // Compute the n-dim box volume in Mpc
    double V =
        std::accumulate(Ls_pc.begin(), Ls_pc.end(), 1, std::multiplies<>());
    V /= 1e6;

    auto sigma = [&](double k) {
        return sqrt(D_init_0 * D_init_0 * V * P(k) / 2);
    };

    deltas = map(deltas, k_mag, [&](std::complex<double> d, double k) {
        // Exact inversion yields transformation formulas for modulus and phase.
        return std::polar(sigma(k) * sqrt(-2 * log(u1())), 2 * M_PI * u2());
    });
    // Note that in the map above we assign Rayleigh + uniform sampled values to
    // all k's. This is wrong for places where k_{x,y,z} = 0, N_{x,y,z}/2
    // (Nyquist frequency). Here, hermitian symmetry enforces the coeffient
    // to be real. Hence we need to sample from a gaussian with twice the
    // standard deviation. Note also that the asignments below is dimension
    // agnostic, i.e. for d=1,2 some coincide.
    auto normal =
        std::bind(std::normal_distribution<>{0, 1}, std::mt19937(rd()));

    for (int i = 0; i < 8; ++i) {
        auto index = [&](int kz, int ky, int kx) {
            return (kz * Ny + ky) * Nx + kx;
        };
        int b1 = (i & 1);
        int b2 = (i & 2) >> 1;
        int b3 = (i & 4) >> 2;
        int idx = index(b3 * Nz / 2, b2 * Ny / 2, b1 * Nx / 2);
        double k = k_mag[idx];
        // rescale unit variance to 2 *sigma(k)
        deltas[idx] = 2 * normal() * sigma(k);
    }

    auto in = reinterpret_cast<fftw_complex*>(deltas.data());
    // Store delta in V for convenience (real vector vs complex vector);
    fftw_plan_ptr c2r(
        fftw_plan_dft_c2r(3, box.Ns.data(), in, state.V.data(), FFTW_ESTIMATE));
    fftw_execute(c2r.get());
    state.V /= V;
}

