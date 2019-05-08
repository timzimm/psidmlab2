#ifndef __CONVOLUTION__
#define __CONVOLUTION__

#include <blaze/math/DynamicVector.h>
#include <fftw3.h>
#include <cassert>
#include <complex>
#include <type_traits>

// Taken from @jeremyfix who took it from GSL. These routines compute the
// optimal padding values to achieve maximum performance of FFTW for DFT
// convolutions.
void factorize(const int n, int* n_factors, int factors[],
               int* implemented_factors);

bool is_optimal(int n, int* implemented_factors);

int find_closest_factor(int n);

// Based on @jeremyfix's convolution benchmark for GSL and FFTW
//
// This workspace is passed around to compute (FFT based) convolutions.
// Typically, such a struct is owned by a class that wishes to perform
// convolutions. Kernels are currently expected to be real.
template <typename T, bool TF = false>
struct convolution_ws {
    int threshold;           // limit up to which the sum approach is used
    bool kernel_up_to_date;  // if false kernel DFT is recomputed
    bool computeLinear;      // if false circular convolution is computed
    bool fast_convolution;   // see first member
    int P;                   // Size of the padded arrays
    int N_offset;  // result sits in signal_padded[N_offset:N_offset + N_signal]
    blaze::DynamicVector<double, TF> kernel_padded;
    blaze::DynamicVector<T, TF> signal_padded;
    blaze::DynamicVector<std::complex<double>, TF>
        kernel_fft;  // holds kernel DFT
    blaze::DynamicVector<std::complex<double>, TF>
        signal_fft;      // holds signal DFT
    fftw_plan r2c;       // r2c plan
    fftw_plan c2r;       // c2r plan
    fftw_plan c2c_for;   // c2c fowards plan
    fftw_plan c2c_back;  // c2c backwards plan

    convolution_ws(const bool computeLinear_, const int N_signal,
                   const int N_kernel);
    ~convolution_ws();
};

template <typename T, bool TF>
convolution_ws<T, TF>::convolution_ws(const bool computeLinear_,
                                      const int N_signal, const int N_kernel)
    : threshold(50),
      kernel_up_to_date(false),
      computeLinear(computeLinear_),
      fast_convolution(N_kernel >= threshold),
      P(0),
      N_offset(computeLinear ? N_kernel / 2 : 0),
      kernel_padded(0),
      signal_padded(0),
      kernel_fft(0),
      signal_fft(0),
      r2c(nullptr),
      c2r(nullptr),
      c2c_for(nullptr),
      c2c_back(nullptr) {
    // DFT based
    if (fast_convolution) {
        P = computeLinear ? find_closest_factor(N_signal + N_kernel / 2)
                          : N_signal;
        // This array holds the convolution result in
        // [N_offset, N_offset + N_signal]
        signal_padded.resize(P);

        auto in = reinterpret_cast<fftw_complex*>(signal_fft.data());

        // Real signal
        if constexpr (std::is_floating_point_v<T>) {
            signal_fft.resize(P / 2 + 1);
            auto out = signal_padded.data();
            // Real signal backwards
            c2r = fftw_plan_dft_c2r_1d(P, in, out, FFTW_ESTIMATE);

            // Complex signal
        } else {
            signal_fft.resize(P);
            auto out = reinterpret_cast<fftw_complex*>(signal_padded.data());
            // Complex signal forwards/backwards
            c2c_for = fftw_plan_dft_1d(P, out, in, FFTW_FORWARD, FFTW_ESTIMATE);
            c2c_back =
                fftw_plan_dft_1d(P, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        }

        // real kernel to complex DFT
        kernel_padded.resize(P);
        kernel_fft.resize(P / 2 + 1);
        r2c = fftw_plan_dft_r2c_1d(
            P, kernel_padded.data(),
            reinterpret_cast<fftw_complex*>(kernel_fft.data()), FFTW_ESTIMATE);

    }
    // Sum based
    else {
        P = computeLinear ? N_signal + N_kernel - 1 : N_signal;
        // This array holds the convolution result in
        // [N_offset, N_offset + N_signal]
        signal_padded.resize(P);
    }
}

template <typename T, bool TF>
convolution_ws<T, TF>::~convolution_ws() {
    fftw_destroy_plan(c2r);
    fftw_destroy_plan(r2c);
    fftw_destroy_plan(c2c_for);
    fftw_destroy_plan(c2c_back);
}

// Applies discrete (linear/circular) convolution by calling one of the routines
// below depending on the threshold value. If N_kernel > threshold DFT based
// fast convolution is used (complexity O(PlogP)). Otherwise the direct
// summation algorithm is used (complexity O(N_signal*N_kernel)). Note that the
// latter is faster for small kernels.
template <typename T, bool TF>
void discrete_convolution(convolution_ws<T, TF>& data,
                          const blaze::DynamicVector<T, TF>& kernel,
                          const blaze::DynamicVector<T, TF>& signal) {
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
        if (data.computeLinear) {
            linear_convolution_sum(data, kernel, signal);
        } else {
            circular_convolution_sum(data, kernel, signal);
        }
    }
}

// Applies sum based linear convolution of signal with kernel
template <typename T, bool TF>
void linear_convolution_sum(convolution_ws<T, TF>& data,
                            const blaze::DynamicVector<T, TF>& kernel,
                            const blaze::DynamicVector<T, TF>& signal) {
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

// Applies sum based circular convolution of signal with kernel.
template <typename T, bool TF>
void circular_convolution_sum(convolution_ws<T, TF>& data,
                              const blaze::DynamicVector<T, TF>& kernel,
                              const blaze::DynamicVector<T, TF>& signal) {
    const int N_kernel = kernel.size();
    const int N_signal = signal.size();

    // Perform convolution of size N_signal
    for (int n = 0; n < data.P; ++n) {
        data.signal_padded[n] = 0;
        for (int k = 0; k < N_kernel; ++k) {
            int i = (n - k) % N_signal;
            i = (i < 0) ? i + N_signal : i;
            data.signal_padded[n] += kernel[k] * signal[i];
        }
    }
}

// Applies DFT based (circular/linear) convolution of signal with kernel
template <typename T, bool TF>
void fast_convolution(convolution_ws<T, TF>& data,
                      const blaze::DynamicVector<double, TF>& kernel,
                      const blaze::DynamicVector<double, TF>& signal) {
    // If the caller changed this switch because kernel changed between calls
    // we need to compute its DFT again after proper zero padding
    if (!data.kernel_up_to_date) {
        // kernel zero padding - bring new kernel into workspace
        blaze::subvector(data.kernel_padded, 0, kernel.size()) = kernel;
        blaze::subvector(data.kernel_padded, kernel.size(),
                         data.P - kernel.size()) = 0;
        fftw_execute(data.r2c);
        data.kernel_up_to_date = true;
    }

    // signal zero padding - bring signal into workspace
    blaze::subvector(data.signal_padded, 0, signal.size()) = signal;
    blaze::subvector(data.signal_padded, signal.size(),
                     data.P - signal.size()) = 0;

    // DFT type depends on type T
    auto sig_fft = reinterpret_cast<fftw_complex*>(data.signal_fft.data());

    if constexpr (std::is_floating_point_v<T>) {
        auto sig = data.signal_padded.data();
        fftw_execute_dft_r2c(data.r2c, sig, sig_fft);
    } else {
        auto sig = reinterpret_cast<fftw_complex*>(data.signal_padded.data());
        fftw_execute_dft(data.c2c_for, sig, sig_fft);
    }

    // Discrete convolution theorem + normalization for IDFT
    data.signal_fft *= 1.0 / data.P * data.kernel_fft;

    if constexpr (std::is_floating_point_v<T>) {
        auto sig = data.signal_padded.data();
        fftw_execute_dft_c2r(data.c2r, sig_fft, sig);

        // TODO Why is this causing a segfault????
        /* fftw_execute(data.c2r); */
    } else {
        auto sig = reinterpret_cast<fftw_complex*>(data.signal_padded.data());
        fftw_execute_dft(data.c2c_back, sig_fft, sig);
    }
}

#endif
