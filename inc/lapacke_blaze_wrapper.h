#ifndef __LAPACKE_WRAP__
#define __LAPACKE_WRAP__
#include <cmath>
#include <complex>
#include "blaze/math/DynamicMatrix.h"
#include "blaze/math/DynamicVector.h"

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

// Forward declare LAPACKE routines
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

// Overloaded wrappers
// TODO if constexpr + is_complex type trait
static void gttrs(int matrix_layout, char trans, int n, int nrhs,
                  const std::complex<double>* dl, const std::complex<double>* d,
                  const std::complex<double>* du,
                  const std::complex<double>* du2, const int* ipiv,
                  std::complex<double>* b, int ldb) {
    LAPACKE_zgttrs(matrix_layout, trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

static void gttrs(int matrix_layout, char trans, int n, int nrhs,
                  const double* dl, const double* d, const double* du,
                  const double* du2, const int* ipiv, double* b, int ldb) {
    LAPACKE_dgttrs(matrix_layout, trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
}

static void gttrf(int n, std::complex<double>* dl, std::complex<double>* d,
                  std::complex<double>* du, std::complex<double>* du2,
                  int* ipiv) {
    LAPACKE_zgttrf(n, dl, d, du, du2, ipiv);
}

static void gttrf(int n, double* dl, double* d, double* du, double* du2,
                  int* ipiv) {
    LAPACKE_dgttrf(n, dl, d, du, du2, ipiv);
}

// Factorizes a general tridiagonal CYCLIC matrix via LU
template <typename T>
static void gctrf(blaze::DynamicVector<T>& dl, blaze::DynamicVector<T>& d,
                  blaze::DynamicVector<T>& du, blaze::DynamicVector<T>& du2,
                  blaze::DynamicVector<int>& ipiv) {
    int N = d.size();
    gttrf(N - 1, dl.data() + 1, d.data(), du.data(), du2.data(), ipiv.data());
}

// Solves a general tridiagonal CYCLIC matrix via the decomposition calculated
// by gctrf
template <typename T>
static void gctrs(const blaze::DynamicVector<T>& dl,
                  const blaze::DynamicVector<T>& d,
                  const blaze::DynamicVector<T>& du,
                  const blaze::DynamicVector<T>& du2,
                  const blaze::DynamicVector<int>& ipiv,
                  blaze::DynamicMatrix<T, blaze::columnMajor>& rhs) {
    using namespace blaze;

    // Size of the cyclic problem
    unsigned int N = rhs.rows();
    unsigned int nrhs = rhs.columns();

    // Condense cyclic matrix into two tridiagonal problems by dropping the last
    // row of the matrix and bringing the last col to the right. Thus, we have
    // an extended right hand side.
    DynamicMatrix<T, columnMajor> rhs_ext(N - 1, nrhs + 1);
    auto x_1 = submatrix(rhs_ext, 0, 0, N - 1, nrhs);
    auto x_2 = column(rhs_ext, nrhs) = 0;

    x_1 = submatrix(rhs, 0, 0, N - 1, nrhs);
    x_2 = -1.0 * dl[0];
    x_2 = -1.0 * du[N - 2];

    // Solve tridiagonal systems.
    // LDB parameter is given by spacing() because blaze adds padding to allow
    // for vectorization
    gttrs(LAPACK_COL_MAJOR, 'N', N - 1, nrhs + 1, dl.data() + 1, d.data(),
          du.data(), du2.data(), ipiv.data(), rhs_ext.data(),
          rhs_ext.spacing());

    // Assemble solution to original problem by linear combination
    auto x_11 = row(x_1, 0);
    auto x_1N = row(x_1, N - 2);

    auto x_s_c = submatrix(rhs, 0, 0, N - 1, nrhs);
    auto x_s_N = row(rhs, N - 1);

    x_s_N -= (du[N - 1] * x_11 + dl[N - 1] * x_1N);
    T normalization =
        1.0 / (d[N - 1] + du[N - 1] * x_2[0] + dl[N - 1] * x_2[N - 2]);
    if (!std::isinf(std::abs(normalization))) x_s_N *= normalization;

    x_s_c = x_1 + x_2 * x_s_N;
}

// Expression templates to enable addition/subtraction of scalars and blaze
// vectors of any type and storage order

template <typename ST>
struct AddScalar {
   public:
    explicit inline AddScalar(ST scalar) : scalar_(scalar) {}

    template <typename T>
    BLAZE_ALWAYS_INLINE decltype(auto) operator()(const T& a) const {
        return a + scalar_;
    }

    template <typename T>
    static constexpr bool simdEnabled() {
        return blaze::HasSIMDAdd<T, ST>::value;
    }

    template <typename T>
    BLAZE_ALWAYS_INLINE decltype(auto) load(const T& a) const {
        BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK(T);
        return a + set(scalar_);
    }

   private:
    ST scalar_;
};
template <typename VT, bool TF, typename Scalar,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<Scalar> > >
decltype(auto) operator+(const blaze::DenseVector<VT, TF>& vec, Scalar scalar) {
    return forEach(~vec, AddScalar<Scalar>(scalar));
}

template <typename Scalar, typename VT, bool TF,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<Scalar> > >
decltype(auto) operator+(Scalar scalar, const blaze::DenseVector<VT, TF>& vec) {
    return forEach(~vec, AddScalar<Scalar>(scalar));
}

template <typename VT, bool TF, typename Scalar,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<Scalar> > >
VT& operator+=(blaze::DenseVector<VT, TF>& vec, Scalar scalar) {
    (~vec) = (~vec) + scalar;
    return ~vec;
}
template <typename ST>
struct SubScalar {
   public:
    explicit inline SubScalar(ST scalar) : scalar_(scalar) {}

    template <typename T>
    BLAZE_ALWAYS_INLINE decltype(auto) operator()(const T& a) const {
        return a - scalar_;
    }

    template <typename T>
    static constexpr bool simdEnabled() {
        return blaze::HasSIMDSub<T, ST>::value;
    }

    template <typename T>
    BLAZE_ALWAYS_INLINE decltype(auto) load(const T& a) const {
        BLAZE_CONSTRAINT_MUST_BE_SIMD_PACK(T);
        return a - set(scalar_);
    }

   private:
    ST scalar_;
};
template <typename VT, bool TF, typename Scalar,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<Scalar> > >
decltype(auto) operator-(const blaze::DenseVector<VT, TF>& vec, Scalar scalar) {
    return forEach(~vec, SubScalar<Scalar>(scalar));
}

template <typename VT, bool TF, typename Scalar,
          typename = blaze::EnableIf_t<blaze::IsNumeric_v<Scalar> > >
VT& operator-=(blaze::DenseVector<VT, TF>& vec, Scalar scalar) {
    (~vec) = (~vec) - scalar;
    return ~vec;
}
#endif
