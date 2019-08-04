#ifndef __LAPACKE_WRAP__
#define __LAPACKE_WRAP__
#include <cmath>
#include <complex>
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

extern "C" {
// ?gttrs solves a tridiagonal matrix system previously factorized with ?gttrf
int LAPACKE_dgttrs(int matrix_layout, char trans, int n, int nrhs,
                   const double* dl, const double* d, const double* du,
                   const double* du2, const int* ipiv, double* b, int ldb);

int LAPACKE_zgttrs(int matrix_layout, char trans, int n, int nrhs,
                   const std::complex<double>* dl,
                   const std::complex<double>* d,
                   const std::complex<double>* du,
                   const std::complex<double>* du2, const int* ipiv,
                   std::complex<double>* b, int ldb);

// ?gttrf factorizes a general tridiagonal matrix via LU
int LAPACKE_dgttrf(int n, double* dl, double* d, double* du, double* du2,
                   int* ipiv);

int LAPACKE_zgttrf(int n, std::complex<double>* dl, std::complex<double>* d,
                   std::complex<double>* du, std::complex<double>* du2,
                   int* ipiv);
}

template <typename T, typename = blaze::EnableIf_t<blaze::IsNumeric_v<T>>>
inline void gttrs(int matrix_layout, char trans, int n, int nrhs, const T* dl,
                  const T* d, const T* du, const T* du2, const int* ipiv, T* b,
                  int ldb) {
    if constexpr (blaze::IsComplexDouble_v<T>)
        LAPACKE_zgttrs(matrix_layout, trans, n, nrhs, dl, d, du, du2, ipiv, b,
                       ldb);
    if constexpr (blaze::IsDouble_v<T>)
        LAPACKE_dgttrs(matrix_layout, trans, n, nrhs, dl, d, du, du2, ipiv, b,
                       ldb);
}

template <typename T, typename = blaze::EnableIf_t<blaze::IsNumeric_v<T>>>
inline void gttrf(int n, T* dl, T* d, T* du, T* du2, int* ipiv) {
    if constexpr (blaze::IsComplexDouble_v<T>)
        LAPACKE_zgttrf(n, dl, d, du, du2, ipiv);
    if constexpr (blaze::IsDouble_v<T>) LAPACKE_dgttrf(n, dl, d, du, du2, ipiv);
}

// Factorizes a general tridiagonal CYCLIC matrix via LU
template <typename T, bool TF,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<T>>>
inline void gctrf(blaze::DynamicVector<T, TF>& dl,
                  blaze::DynamicVector<T, TF>& d,
                  blaze::DynamicVector<T, TF>& du,
                  blaze::DynamicVector<T, TF>& du2,
                  blaze::DynamicVector<int, TF>& ipiv) {
    int N = d.size();
    gttrf(N - 1, dl.data() + 1, d.data(), du.data(), du2.data(), ipiv.data());
}

// Solves a general tridiagonal CYCLIC matrix via the decomposition calculated
// by gctrf
template <typename T, bool TF, bool SO,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<T>>>
inline void gctrs(const blaze::DynamicVector<T, TF>& dl,
                  const blaze::DynamicVector<T, TF>& d,
                  const blaze::DynamicVector<T, TF>& du,
                  const blaze::DynamicVector<T, TF>& du2,
                  const blaze::DynamicVector<int, TF>& ipiv,
                  blaze::DynamicMatrix<T, SO>& rhs) {
    using namespace blaze;

    const int matrix_layout =
        (SO == columnMajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
    // Size of the cyclic problem
    const int N = rhs.rows();
    const int nrhs = rhs.columns();

    // Condense cyclic matrix into two tridiagonal problems by dropping the last
    // row of the matrix and bringing the last col to the right. Thus, we have
    // an extended right hand side.
    rhs.resize(N, nrhs + 1);

    auto x_2 = subvector(column(rhs, nrhs), 0, N - 1) = T{};
    x_2[0] = -1.0 * dl[0];
    x_2[N - 2] = -1.0 * du[N - 2];

    // Solve tridiagonal systems by calling LAPACKE
    // LDB parameter is given by spacing() because blaze adds padding to allow
    // for vectorization
    gttrs(matrix_layout, 'N', N - 1, nrhs + 1, dl.data() + 1, d.data(),
          du.data(), du2.data(), ipiv.data(), rhs.data(), rhs.spacing());

    // Assemble solution to original problem by linear combination
    auto x_1 = submatrix(rhs, 0, 0, N - 1, nrhs);
    auto x_11 = row(x_1, 0);
    auto x_1N = row(x_1, N - 2);

    auto x_s_N = subvector(row(rhs, N - 1), 0, nrhs);

    x_s_N -= (du[N - 1] * x_11 + dl[N - 1] * x_1N);
    T normalization =
        1.0 / (d[N - 1] + du[N - 1] * x_2[0] + dl[N - 1] * x_2[N - 2]);
    if (!std::isinf(std::abs(normalization))) x_s_N *= normalization;

    x_1 += x_2 * x_s_N;
    rhs.resize(N, nrhs);
}

#endif
