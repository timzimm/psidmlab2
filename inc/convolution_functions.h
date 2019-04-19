#ifndef __CONVOLUTION__
#define __CONVOLUTION__

#include <blaze/math/DynamicVector.h>
#include <fftw3.h>
#include <cassert>

// Based on @jeremyfix's convolution benchmark for GSL and FFTW
//
// Functions performing real data convolutions. It is somewhat odd to wrap this
// into a class. Hence, we employ a more functional/ C-like approach until I get
// sufficiently annoyed about this.

// This workspace is passed around to compute (FFT based) convolutions.
// Typically, such a struct is owned by a class that wishes to perform
// convolutions.
struct convolution_ws {
    int threshold;           // limit up to which the sum approach is used
    bool kernel_up_to_date;  // if false kernel DFT is recomputed
    bool computeLinear;      // if false circular convolution is computed
    bool fast_convolution;   // see first member
    int P;                   // Size of the padded arrays
    int N_offset;  // result sits in signal_padded[N_offset:N_offset + N_signal]
    blaze::DynamicVector<double> kernel_padded;
    blaze::DynamicVector<double> signal_padded;
    blaze::DynamicVector<std::complex<double>> kernel_fft;  // holds kernel DFT
    blaze::DynamicVector<std::complex<double>> signal_fft;  // holds signal DFT
    fftw_plan forward;                                      // r2c plan
    fftw_plan backward;                                     // c2r plan

    convolution_ws(const bool computeLinear_, const int N_signal,
                   const int N_kernel);
    ~convolution_ws();
};

// Applies discrete (linear/circular) convolution by calling one of the routines
// below depending on the threshold value. If N_kernel > threshold DFT based
// fast convolution is used (complexity O(PlogP)). Otherwise the direct
// summation algorithm is used (complexity O(N_signal*N_kernel)). Note that the
// latter is faster for small kernels.
void discrete_convolution(convolution_ws& data,
                          const blaze::DynamicVector<double>& kernel,
                          const blaze::DynamicVector<double>& signal);

// Applies DFT convolution (linear/circular) on signal with kernel using
// allocated memeory in ws. If data.computeLinear is false, the shift theorem is
// applied to recenter the output data. This might change in the future.
void fast_convolution(convolution_ws& data,
                      const blaze::DynamicVector<double>& kernel,
                      const blaze::DynamicVector<double>& signal);

// Applies sum based linear convolution of signal with kernel and recenters its
// output data.signal_padded on the fly.
void linear_convolution_sum(convolution_ws& data,
                            const blaze::DynamicVector<double>& kernel,
                            const blaze::DynamicVector<double>& signal);

// Applies sum based circular convolution of signal with kernel.
void circular_convolution_sum(convolution_ws& data,
                              const blaze::DynamicVector<double>& kernel,
                              const blaze::DynamicVector<double>& signal);

// Taken from @jeremyfix who took it from GSL. These routines compute the
// optimal padding values to achieve maximum performance of FFTW for DFT
// convolutions.
void factorize(const int n, int* n_factors, int factors[],
               int* implemented_factors);

bool is_optimal(int n, int* implemented_factors);

int find_closest_factor(int n);

#endif
