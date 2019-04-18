#include "convolution_functions.h"
#include <iostream>

convolution_ws::convolution_ws(const bool computeLinear_, const int N_signal,
                               const int N_kernel)
    : threshold(20),
      kernel_up_to_date(false),
      computeLinear(computeLinear_),
      fast_convolution(N_kernel >= threshold),
      P(0),
      N_offset(computeLinear ? N_kernel / 2 : 0),
      kernel_padded(0),
      signal_padded(0),
      kernel_fft(0),
      signal_fft(0),
      forward(nullptr),
      backward(nullptr) {
    if (fast_convolution) {
        P = computeLinear ? find_closest_factor(N_signal + N_kernel / 2)
                          : N_signal;
        // We only allocate extra memory for DFTs if we really perform them
        signal_fft.resize(P / 2 + 1);
        kernel_fft.resize(P / 2 + 1);
        kernel_padded.resize(P);

        // Domain and target vectors might change. Plans, however, are still
        // valid in that case
        forward = fftw_plan_dft_r2c_1d(
            P, kernel_padded.data(),
            reinterpret_cast<fftw_complex*>(kernel_fft.data()), FFTW_ESTIMATE);
        backward = fftw_plan_dft_c2r_1d(
            P, reinterpret_cast<fftw_complex*>(kernel_fft.data()),
            kernel_padded.data(), FFTW_ESTIMATE);
    } else {
        P = computeLinear ? N_signal + N_kernel - 1 : N_signal;
    }
    // This array holds the convolution result in
    // [N_offset, N_offset + N_signal]
    signal_padded.resize(P);
}

convolution_ws::~convolution_ws() {
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
}

void discrete_convolution(convolution_ws& data,
                          const blaze::DynamicVector<double>& kernel,
                          const blaze::DynamicVector<double>& signal) {
    const int N_kernel = kernel.size();
    const int N_signal = signal.size();

    // Check if input sizes make sense
    assert(N_kernel > 0 || N_signal > 0 || N_signal >= N_kernel);

    // Dispatch correct convolution routine (i.e. sum or DFT based)
    if (data.fast_convolution) {
        // DFT based convolutions handle both linear and circular cases by
        // appropriate padding P
        fast_convolution(data, kernel, signal);
    } else {
        if (data.computeLinear)
            linear_convolution_sum(data, kernel, signal);
        else
            circular_convolution_sum(data, kernel, signal);
    }
}
void fast_convolution(convolution_ws& data,
                      const blaze::DynamicVector<double>& kernel,
                      const blaze::DynamicVector<double>& signal) {
    // If the caller changed this switch because kernel changed between calls
    // we need to compute its DFT again after proper zero padding
    if (!data.kernel_up_to_date) {
        // kernel zero padding
        blaze::subvector(data.kernel_padded, 0, kernel.size()) = kernel;
        blaze::subvector(data.kernel_padded, kernel.size(),
                         data.P - kernel.size()) = 0;
        fftw_execute_dft_r2c(
            data.forward, data.kernel_padded.data(),
            reinterpret_cast<fftw_complex*>(data.kernel_fft.data()));
        data.kernel_up_to_date = true;
    }

    // signal zero padding
    blaze::subvector(data.signal_padded, 0, signal.size()) = signal;
    blaze::subvector(data.signal_padded, signal.size(),
                     data.P - signal.size()) = 0;

    // Take DFT of input signal
    fftw_execute_dft_r2c(
        data.forward, data.signal_padded.data(),
        reinterpret_cast<fftw_complex*>(data.signal_fft.data()));

    // Discrete convolution theorem + normalization for IDFT
    data.signal_fft *= 1.0 / data.P * data.kernel_fft;

    fftw_execute_dft_r2c(
        data.backward, data.signal_padded.data(),
        reinterpret_cast<fftw_complex*>(data.signal_fft.data()));
}

void linear_convolution_sum(convolution_ws& data,
                            const blaze::DynamicVector<double>& kernel,
                            const blaze::DynamicVector<double>& signal) {
    const int N_kernel = kernel.size();
    const int N_signal = signal.size();

    // Perform FULL convolution of size P = N_signal + N_kernel - 1
    for (int n = 0; n < data.P; ++n) {
        data.signal_padded[n] = 0;
        int khigh = std::min(N_signal - 1, n);
        int klow = std::max(0, n - N_kernel + 1);

        for (int k = klow; k <= khigh; ++k) {
            data.signal_padded[n] += signal[k] * kernel[n - k];
        }
    }
}

void circular_convolution_sum(convolution_ws& data,
                              const blaze::DynamicVector<double>& kernel,
                              const blaze::DynamicVector<double>& signal) {
    const int N_kernel = kernel.size();
    const int N_signal = signal.size();

    // Perform convolution of size N_signal
    for (int n = 0; n < data.P; ++n) {
        for (int k = 0; k < N_kernel; ++k) {
            int i = (n - k) % N_signal;
            i = (i < 0) ? i + N_signal : i;
            data.signal_padded[n] += kernel[k] * signal[i];
        }
    }
}

void factorize(const int n, int* n_factors, int factors[],
               int* implemented_factors) {
    int nf = 0;
    int ntest = n;
    int factor;
    int i = 0;

    if (n == 1) {
        factors[0] = 1;
        *n_factors = 1;
        return;
    }

    while (implemented_factors[i] && ntest != 1) {
        factor = implemented_factors[i];
        while ((ntest % factor) == 0) {
            ntest = ntest / factor;
            factors[nf] = factor;
            nf++;
        }
        i++;
    }

    if (ntest != 1) {
        factors[nf] = ntest;
        nf++;
    }

    {
        int product = 1;

        for (i = 0; i < nf; i++) {
            product *= factors[i];
        }
    }

    *n_factors = nf;
}

bool is_optimal(int n, int* implemented_factors) {
    // We check that n is not a multiple of 4*4*4*2
    if (n % 4 * 4 * 4 * 2 == 0) return false;

    int nf = 0;
    int factors[64];
    int i = 0;
    factorize(n, &nf, factors, implemented_factors);

    while (implemented_factors[i]) {
        if (factors[nf - 1] == implemented_factors[i]) return true;
        ++i;
    }
    return false;
}

int find_closest_factor(int n) {
    // FFTW optimizes for these factors
    int implemented_factors[] = {13, 11, 7, 5, 3, 2, 0};
    int j;
    if (is_optimal(n, implemented_factors)) return n;
    j = n + 1;
    while (!is_optimal(j, implemented_factors)) ++j;
    return j;
}
