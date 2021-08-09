#pragma once
#pragma warning (disable:4819)

//#pragma comment(lib, "mkl_rt.lib")


#include "mkl.h"
#include <omp.h>


class BLASInterface
{
private:
    static int s_NumThread;

public:
    static void SetNumThreads(int nt);
    static inline int GetNumThreads() { return s_NumThread; }

    // [x], [y]: Vector x, y
    // [A], [B], [C]: matrix A, B, C
    // alpha, beta: scalar

    // BLAS Level 1
    // Dot Product: result = [x] * [y]
    template<typename DataType>
    static DataType DOT(const int n, const DataType* x, const int incx, const DataType* y, const int incy);

    // BLAS Level 1
    // Vector-Scalar Product: [x] = alpha * [x]
    template<typename DataType>
    static void SCAL(const int n, const DataType alpha, DataType* x, const int incx);

    // BLAS Level 1
    // Vector-Scalar Product + Vector: [y] = alpha * [x] + [y]
    template<typename DataType>
    static void AXPY(const int n, const DataType alpha, 
        const DataType* x, const int incx, DataType* y, const int incy);

    // BLAS Level 2
    // [y] = alpha * [A] * [x] + beta * [y];  transA = CBlasNoTrans
    // [y] = alpha * [A]' * [x] + beta * [y];  transA = CBlasTrans
    // [y] = alpha * conjugate([A]') * [x] + beta * [y];  transA = CBlasConjTrans
    // [A]: m by n
    // CBLAS_LAYOUT { CblasRowMajor, CblasColMajor }
    // CBLAS_TRANSPOSE { CblasNoTrans, CblasTrans, CblasConjTrans }
    // lda: Leading dimension of [A]. Row-Major n, Column-Major m
    // k: Row-Major m, Column-Major n
    template<typename DataType>
    static void GEMV(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const int m, const int n,
        const DataType alpha, const DataType* A, const int lda, const DataType* x, const int incx,
        const DataType beta, DataType* y, const int incy);
    
    // BLAS Level 3
    // [C] = alpha * op([A]) * op([B]) + beta * [C]
    // op([A]): m by k. [A], [A]', conjugate([A]')
    // op([B]): k by n. [B], [B]', conjugate([B]')
    // [C]: m by n.
    // lda: Leading dimension of [A].
    // ldb: Leading dimension of [B].
    // ldc: Leading demension of [C].
    template<typename DataType>
    static void GEMM(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, 
        const int m, const int n, const int k, const DataType alpha, const DataType* A, const int lda, 
        const DataType* B, const int ldb, const DataType beta, DataType* C, const int ldc);

    // BLAS Level 3
    // op([A]) * [X] = alpha * [B]: CBLAS_SIDE = CblasLeft
    // [X] * op([A]) = alpha * [B]: CBLAS_SIDE = CblasRight
    // CBLAS_UPLO { CblasUpper, CblasLower }
    // [A]: Unit or non-unit, upper or lower triangular matrix
    // [X]: Unknown matrix. m by n
    // [B]: m by n
    // CBLAS_DIAG { CblasNonUnit, CblasUnit }
    // Unit triangular matrix: Triangular matrix with 1 on the diagonal
    template<typename DataType>
    static void TRSM(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, 
        const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int m, const int n, 
        const DataType alpha, const DataType* A, const int lda, DataType* B, const int ldb);

    // Linear Equations LAPACK
    // Computes the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix
    // [A] = [U]' * [U]: uplo = 'U'
    // [A] = [L] * [L]': uplo = 'L'
    // matrix_layout { LAPACK_ROW_MAJOR, LAPACK_COL_MAJOR }
    // uplo { 'U', 'L' }
    // Supports Progress Function (mkl_progress)
    // Important Node: Result are overwritten in [A].
    // For upper case, lower parts should be filled with zero after calculation.
    // For lower case, upper parts should be filled with zero after calculation. 
    template<typename DataType>
    static int POTRF(int matrix_layout, char uplo, int n, DataType* A, int lda);

    // Linear Equations LAPACK
    // Solves a system of linear equations with a Cholesky-factored 
    // symmetric (Hermitian) positive-definite coefficient matrix
    // Solve [A][x] = [B]
    // nrhs: The number of right-hand sides. If [B] is a column vector, this value is 1.
    // Important Node: POTRF should be called before calling this function.
    template<typename DataType>
    static int POTRS(int matrix_layout, char uplo, int n, int nrhs, 
        const DataType* A, int lda, DataType* B, int ldb);

    // Linear Equations LAPACK
    // Computes the Bunch-Kaufman factorization of a symmetric matrix
    // [A] = [U] * [D] * [U]': uplo = 'U'
    // [A] = [L] * [D] * [L]': uplo = 'L'
    // Supports Progress Function (mkl_progress)
    // ipiv: Array at least n-size. Contains details of the interchanges and the block structure of D.
    template<typename DataType>
    static int SYTRF(int matrix_layout, char uplo, int n, DataType* A, int lda, int* ipiv);
};


inline void BLASInterface::SetNumThreads(int nt)
{
    s_NumThread = nt;
    MKL_Set_Num_Threads(nt);
    omp_set_num_threads(nt);
}

template<>
inline float BLASInterface::DOT(const int n, const float* x, const int incx, const float* y, const int incy)
{
    return cblas_sdot(n, x, incx, y, incy);
}

template<>
inline double BLASInterface::DOT(const int n, const double* x, const int incx, const double* y, const int incy)
{
    return cblas_ddot(n, x, incx, y, incy);
}

template<typename DataType>
inline DataType BLASInterface::DOT(const int n, const DataType* x, const int incx, const DataType* y, const int incy)
{
    assert(false);
    return static_cast<DataType>(0.0);
}

template<>
inline void BLASInterface::SCAL(const int n, const float alpha, float* x, const int incx)
{
    cblas_sscal(n, alpha, x, incx);
}

template<>
inline void BLASInterface::SCAL(const int n, const double alpha, double* x, const int incx)
{
    cblas_dscal(n, alpha, x, incx);
}

template<typename DataType>
inline void BLASInterface::SCAL(const int n, const DataType alpha, DataType* x, const int incx)
{
    assert(false);
}

template<>
inline void BLASInterface::AXPY(const int n, const float alpha, 
    const float* x, const int incx, float* y, const int incy)
{
    cblas_saxpy(n, alpha, x, incx, y, incy);
}

template<>
inline void BLASInterface::AXPY(const int n, const double alpha, 
    const double* x, const int incx, double* y, const int incy)
{
    cblas_daxpy(n, alpha, x, incx, y, incy);
}

template<typename DataType>
inline void BLASInterface::AXPY(const int n, const DataType alpha, 
    const DataType* x, const int incx, DataType* y, const int incy)
{
    assert(false);
}

template<>
inline void BLASInterface::GEMV(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const int m, const int n,
    const float alpha, const float* A, const int lda, const float* x, const int incx,
    const float beta, float* y, const int incy)
{
    cblas_sgemv(layout, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template<>
inline void BLASInterface::GEMV(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const int m, const int n,
    const double alpha, const double* A, const int lda, const double* x, const int incx,
    const double beta, double* y, const int incy)
{
    cblas_dgemv(layout, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template<typename DataType>
inline void BLASInterface::GEMV(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const int m, const int n,
    const DataType alpha, const DataType* A, const int lda, const DataType* x, const int incx,
    const DataType beta, DataType* y, const int incy)
{
    assert(false);
}

template<>
inline void BLASInterface::GEMM(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
    const int m, const int n, const int k, const float alpha, const float* A, const int lda,
    const float* B, const int ldb, const float beta, float* C, const int ldc)
{
    cblas_sgemm(layout, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<>
inline void BLASInterface::GEMM(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
    const int m, const int n, const int k, const double alpha, const double* A, const int lda,
    const double* B, const int ldb, const double beta, double* C, const int ldc)
{
    cblas_dgemm(layout, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

template<typename DataType>
inline void BLASInterface::GEMM(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
    const int m, const int n, const int k, const DataType alpha, const DataType* A, const int lda,
    const DataType* B, const int ldb, const DataType beta, DataType* C, const int ldc)
{
    assert(false);
}

template<>
inline void BLASInterface::TRSM(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int m, const int n,
    const float alpha, const float* A, const int lda, float* B, const int ldb)
{
    cblas_strsm(layout, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

template<>
inline void BLASInterface::TRSM(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int m, const int n,
    const double alpha, const double* A, const int lda, double* B, const int ldb)
{
    cblas_dtrsm(layout, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

template<typename DataType>
inline void BLASInterface::TRSM(const CBLAS_LAYOUT layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE transA, const CBLAS_DIAG diag, const int m, const int n,
    const DataType alpha, const DataType* A, const int lda, DataType* B, const int ldb)
{
    assert(false);
}

template<>
inline int BLASInterface::POTRF(int matrix_layout, char uplo, int n, float* A, int lda)
{
    return LAPACKE_spotrf(matrix_layout, uplo, n, A, lda);
}

template<>
inline int BLASInterface::POTRF(int matrix_layout, char uplo, int n, double* A, int lda)
{
    return LAPACKE_dpotrf(matrix_layout, uplo, n, A, lda);
}

template<typename DataType>
inline int BLASInterface::POTRF(int matrix_layout, char uplo, int n, DataType* A, int lda)
{
    assert(false);
    return -1;
}

template<>
inline int BLASInterface::POTRS(int matrix_layout, char uplo, int n, int nrhs,
    const float* A, int lda, float* B, int ldb)
{
    return LAPACKE_spotrs(matrix_layout, uplo, n, nrhs, A, lda, B, ldb);
}

template<>
inline int BLASInterface::POTRS(int matrix_layout, char uplo, int n, int nrhs,
    const double* A, int lda, double* B, int ldb)
{
    return LAPACKE_dpotrs(matrix_layout, uplo, n, nrhs, A, lda, B, ldb);
}

template<typename DataType>
inline int BLASInterface::POTRS(int matrix_layout, char uplo, int n, int nrhs,
    const DataType* A, int lda, DataType* B, int ldb)
{
    assert(false);
    return -1;
}

template<>
inline int BLASInterface::SYTRF(int matrix_layout, char uplo, int n, float* A, int lda, int* ipiv)
{
    return LAPACKE_ssytrf(matrix_layout, uplo, n, A, lda, ipiv);
}

template<>
inline int BLASInterface::SYTRF(int matrix_layout, char uplo, int n, double* A, int lda, int* ipiv)
{
    return LAPACKE_dsytrf(matrix_layout, uplo, n, A, lda, ipiv);
}

template<typename DataType>
inline int BLASInterface::SYTRF(int matrix_layout, char uplo, int n, DataType* A, int lda, int* ipiv)
{
    assert(false);
    return -1;
}
