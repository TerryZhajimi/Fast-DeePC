
#include <math.h>
#include <ctype.h>
#include "deepc_boxqp_optimized.c"

static int is_t(const char *c) {
    return c && (c[0] == 't' || c[0] == 'T' || c[0] == 'c' || c[0] == 'C');
}

void dcopy_(const long *n, const double *x, const long *incx, double *y, const long *incy) {
    long i;
    for (i = 0; i < *n; ++i) y[i * (*incy)] = x[i * (*incx)];
}

void daxpy_(const long *n, const double *alpha, const double *x, const long *incx, double *y, const long *incy) {
    long i;
    for (i = 0; i < *n; ++i) y[i * (*incy)] += (*alpha) * x[i * (*incx)];
}

void dscal_(const long *n, const double *alpha, double *x, const long *incx) {
    long i;
    for (i = 0; i < *n; ++i) x[i * (*incx)] *= *alpha;
}

double ddot_(const long *n, const double *x, const long *incx, const double *y, const long *incy) {
    long i;
    double s = 0.0;
    for (i = 0; i < *n; ++i) s += x[i * (*incx)] * y[i * (*incy)];
    return s;
}

void dgemv_(const char *trans, const long *m, const long *n, const double *alpha,
            const double *A, const long *lda, const double *x, const long *incx,
            const double *beta, double *y, const long *incy) {
    long i, j;
    if (!is_t(trans)) {
        for (i = 0; i < *m; ++i) {
            double s = 0.0;
            for (j = 0; j < *n; ++j) s += A[i + j * (*lda)] * x[j * (*incx)];
            y[i * (*incy)] = (*alpha) * s + (*beta) * y[i * (*incy)];
        }
    } else {
        for (j = 0; j < *n; ++j) {
            double s = 0.0;
            for (i = 0; i < *m; ++i) s += A[i + j * (*lda)] * x[i * (*incx)];
            y[j * (*incy)] = (*alpha) * s + (*beta) * y[j * (*incy)];
        }
    }
}

void dgemm_(const char *transa, const char *transb, const long *m, const long *n, const long *k,
            const double *alpha, const double *A, const long *lda, const double *B, const long *ldb,
            const double *beta, double *C, const long *ldc) {
    long i, j, l;
    int ta = is_t(transa);
    int tb = is_t(transb);
    for (j = 0; j < *n; ++j) {
        for (i = 0; i < *m; ++i) {
            double s = 0.0;
            for (l = 0; l < *k; ++l) {
                double a = ta ? A[l + i * (*lda)] : A[i + l * (*lda)];
                double b = tb ? B[j + l * (*ldb)] : B[l + j * (*ldb)];
                s += a * b;
            }
            C[i + j * (*ldc)] = (*alpha) * s + (*beta) * C[i + j * (*ldc)];
        }
    }
}

void dpotrf_(const char *uplo, const long *n, double *A, const long *lda, long *info) {
    long i, j, k;
    *info = 0;
    if (uplo && (uplo[0] == 'u' || uplo[0] == 'U')) {
        *info = -1;
        return;
    }
    for (j = 0; j < *n; ++j) {
        double diag = A[j + j * (*lda)];
        for (k = 0; k < j; ++k) diag -= A[j + k * (*lda)] * A[j + k * (*lda)];
        if (!(diag > 0.0) || !isfinite(diag)) {
            *info = j + 1;
            return;
        }
        A[j + j * (*lda)] = sqrt(diag);
        for (i = j + 1; i < *n; ++i) {
            double s = A[i + j * (*lda)];
            for (k = 0; k < j; ++k) s -= A[i + k * (*lda)] * A[j + k * (*lda)];
            A[i + j * (*lda)] = s / A[j + j * (*lda)];
        }
        for (k = j + 1; k < *n; ++k) A[j + k * (*lda)] = 0.0;
    }
}

void dpotrs_(const char *uplo, const long *n, const long *nrhs, const double *A, const long *lda,
             double *B, const long *ldb, long *info) {
    long i, j, k;
    *info = 0;
    if (uplo && (uplo[0] == 'u' || uplo[0] == 'U')) {
        *info = -1;
        return;
    }
    for (j = 0; j < *nrhs; ++j) {
        for (i = 0; i < *n; ++i) {
            double s = B[i + j * (*ldb)];
            for (k = 0; k < i; ++k) s -= A[i + k * (*lda)] * B[k + j * (*ldb)];
            B[i + j * (*ldb)] = s / A[i + i * (*lda)];
        }
        for (i = *n - 1; i >= 0; --i) {
            double s = B[i + j * (*ldb)];
            for (k = i + 1; k < *n; ++k) s -= A[k + i * (*lda)] * B[k + j * (*ldb)];
            B[i + j * (*ldb)] = s / A[i + i * (*lda)];
        }
    }
}
