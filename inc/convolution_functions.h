#ifndef __CONVOLUTION__
#define __CONVOLUTION__

#include <blaze/math/DynamicVector.h>
#include <blaze/util/TypeTraits.h>
#include <fftw3.h>
#include <complex>

// Taken from from GSL. These routines compute the
// optimal padding values to achieve maximum performance of FFTW for DFT
// convolutions.
void factorize(const int n, int* n_factors, int factors[],
               int* implemented_factors);

bool is_optimal(int n, int* implemented_factors);

int find_closest_factor(int n);

// This workspace is passed around to compute (FFT based) convolutions.
// Typically, such a struct is owned by a class that wishes to perform
// convolutions. Kernels are currently expected to be real and one-dimensional
template <typename T, bool TF = false>
struct convolution_ws {
    const int N_signal;
    const int N_kernel;
    const int threshold;          // limit up to which the sum approach is used
    bool kernel_up_to_date;       // if false kernel DFT is recomputed
    const bool computeLinear;     // if false circular convolution is computed
    const bool fast_convolution;  // see first member
    int P;                        // Size of the padded arrays
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
                                      const int N_signal_, const int N_kernel_)
    : N_signal(N_signal_),
      N_kernel(N_kernel_),
      threshold(50),
      kernel_up_to_date(false),
      computeLinear(computeLinear_),
      fast_convolution(N_kernel >= threshold),
      P(0),
      kernel_padded(N_kernel),
      signal_padded(N_signal),
      kernel_fft(0),
      signal_fft(0),
      r2c(nullptr),
      c2r(nullptr),
      c2c_for(nullptr),
      c2c_back(nullptr) {
    static_assert(blaze::IsFloatingPoint_v<T> || blaze::IsComplex_v<T>);
    // DFT based
    if (fast_convolution) {
        P = computeLinear ? find_closest_factor(N_signal + N_kernel / 2)
                          : N_signal;

        auto in = reinterpret_cast<fftw_complex*>(signal_fft.data());

        // Real signal
        if constexpr (blaze::IsFloatingPoint_v<T>) {
            signal_fft.reserve(P / 2 + 1);
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
        kernel_fft.reserve(P / 2 + 1);
        kernel_fft.resize(P / 2 + 1);

        r2c = fftw_plan_dft_r2c_1d(
            P, kernel_padded.data(),
            reinterpret_cast<fftw_complex*>(kernel_fft.data()), FFTW_ESTIMATE);

    }
    // Sum based
    else {
        P = computeLinear ? N_signal + N_kernel - 1 : N_signal;
    }
    // Allow for efficient resize
    signal_padded.reserve(P);
    kernel_padded.reserve(P);
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
void discrete_convolution(convolution_ws<T, TF>& data) {
    // Add zero padding
    data.signal_padded.resize(data.P);
    data.kernel_padded.resize(data.P);

    // Dispatch correct convolution routine (i.e. sum or DFT based)
    if (data.fast_convolution) {
        // DFT based convolutions handle both linear and circular cases by
        // appropriate padding P
        fast_convolution(data);
    } else {
        if (data.computeLinear) {
            linear_convolution_sum(data);
        } else {
            circular_convolution_sum(data);
        }
    }

    std::rotate(data.signal_padded.begin(),
                data.signal_padded.begin() + data.N_kernel / 2,
                data.signal_padded.end());

    data.signal_padded.resize(data.N_signal);
    data.kernel_padded.resize(data.N_kernel);
}

// Applies sum based linear convolution of signal with kernel
template <typename T, bool TF>
void linear_convolution_sum(convolution_ws<T, TF>& ws) {
    const blaze::DynamicVector<T, TF> signal =
        blaze::subvector(ws.signal_padded, 0, ws.N_signal);

    // Perform FULL convolution of size P = N_signal + N_kernel - 1
    for (int n = 0; n < ws.P; ++n) {
        ws.signal_padded[n] = 0;
        int khigh = std::min(ws.N_signal - 1, n);
        int klow = std::max(0, n - ws.N_kernel + 1);

        for (int k = klow; k <= khigh; ++k) {
            ws.signal_padded[n] += signal[k] * ws.kernel_padded[n - k];
        }
    }
}

// Applies sum based circular convolution of signal with kernel.
template <typename T, bool TF>
void circular_convolution_sum(convolution_ws<T, TF>& ws) {
    const blaze::DynamicVector<T, TF> signal =
        blaze::subvector(ws.signal_padded, 0, ws.N_signal);

    // Perform convolution of size N_signal
    for (int n = 0; n < ws.P; ++n) {
        ws.signal_padded[n] = 0;
        for (int k = 0; k < ws.N_kernel; ++k) {
            int i = (n - k) % ws.N_signal;
            i = (i < 0) ? i + ws.N_signal : i;
            ws.signal_padded[n] += ws.kernel_padded[k] * signal[i];
        }
    }
}

// Applies DFT based (circular/linear) convolution of signal with kernel
template <typename T, bool TF>
void fast_convolution(convolution_ws<T, TF>& data) {
    // If the caller changed this switch because kernel changed between calls
    // we need to compute its DFT again after proper zero padding
    if (!data.kernel_up_to_date) {
        auto in = data.kernel_padded.data();
        auto out = reinterpret_cast<fftw_complex*>(data.kernel_fft.data());
        fftw_execute_dft_r2c(data.r2c, in, out);
        data.kernel_up_to_date = true;
    }

    // DFT type depends on type T
    auto sig_fft = reinterpret_cast<fftw_complex*>(data.signal_fft.data());

    if constexpr (blaze::IsFloatingPoint_v<T>) {
        auto sig = data.signal_padded.data();
        fftw_execute_dft_r2c(data.r2c, sig, sig_fft);
    } else {
        auto sig = reinterpret_cast<fftw_complex*>(data.signal_padded.data());
        fftw_execute_dft(data.c2c_for, sig, sig_fft);
    }

    if constexpr (blaze::IsFloatingPoint_v<T>) {
        // Discrete convolution theorem + normalization
        data.signal_fft *= 1.0 / data.P * data.kernel_fft;
        auto sig = data.signal_padded.data();
        fftw_execute_dft_c2r(data.c2r, sig_fft, sig);
    } else {
        const int N = data.P;
        // Discrete convolution theorem + hermitian property + normalization
        subvector(data.signal_fft, 0, N / 2 + 1) *= data.kernel_fft;
        subvector(data.signal_fft, N / 2 + 1, N - (N / 2 + 1)) *=
            conj(reverse(subvector(data.kernel_fft, 1, N - (N / 2 + 1))));
        data.signal_fft *= 1.0 / data.P;
        auto sig = reinterpret_cast<fftw_complex*>(data.signal_padded.data());
        fftw_execute_dft(data.c2c_back, sig_fft, sig);
    }
}

#endif
