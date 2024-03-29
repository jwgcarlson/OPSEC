#ifndef OPSEC_H
#define OPSEC_H

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "cfg.h"

#ifdef __cplusplus
extern "C" {
#endif

/* OPSEC globals */
extern int opsec_verbose;

/* Parse command line and fill in a Config object.  This Config object must
 * be cfg_destroy()'ed when no longer needed. */
Config opsec_init(int argc, char * const argv[], const char *usage);

/* Exit OPSEC safely and cleanly. */
void opsec_exit(int status);

/* Abort OPSEC abruptly and loudly. */
void opsec_abort(int status);

/* Convenience routines for logging status and error messages:
 *   opsec_info: general status messages (stdout, always)
 *   opsec_debug: debugging-related messages (stdout, OPSEC_DEBUG only)
 *   opsec_warn: warning messages (stderr, OPSEC_WARN only)
 *   opsec_error: error messages (stderr, always)
 * These routines are to be used as replacements for printf.  For MPI programs,
 * messages will only be printed on the root process. */
int opsec_info(const char *format, ...);
int opsec_debug(const char *format, ...);
int opsec_warn(const char *format, ...);
int opsec_error(const char *format, ...);

/* Generic print function used by the above convenience routines. */
int opsec_vfprintf(FILE* stream, const char *format, va_list ap);


#ifdef OPSEC_DEBUG
    /* Verbose wrappers around standard memory allocation routines */
#   define opsec_malloc(size) opsec_debug_malloc(size, __FILE__, __LINE__)
#   define opsec_calloc(nmemb, size) opsec_debug_calloc(nmemb, size, __FILE__, __LINE__)

    /* Assertions */
#   define opsec_assert(expr) assert(expr)

#else
#   define opsec_malloc(size) malloc(size)
#   define opsec_calloc(nmemb, size) calloc(nmemb, size)
#   define opsec_assert(expr) do {} while(0)
#endif


/* Choose between single and double precision */
#ifdef OPSEC_SINGLE
    typedef float real;
#   define REAL_MPI_TYPE MPI_FLOAT
#   define blas_dot F77_FUNC(sdot, SDOT)
#   define blas_gemv F77_FUNC(sgemv, SGEMV)
#   define blas_gemm F77_FUNC(sgemm, SGEMM)
#   define blas_axpy F77_FUNC(saxpy, SAXPY)
#   define lapack_syev F77_FUNC(ssyev, SSYEV)
#   define lapack_gesvd F77_FUNC(sgesvd, SGESVD)
#   define lapack_getrf F77_FUNC(sgetrf, SGETRF)
#   define lapack_getri F77_FUNC(sgetri, SGETRI)
#   define arpack_saupd F77_FUNC(ssaupd, SSAUPD)
#   define arpack_seupd F77_FUNC(sseupd, SSEUPD)
#   define parpack_psaupd F77_FUNC(pssaupd, PSSAUPD)
#   define parpack_pseupd F77_FUNC(psseupd, PSSEUPD)
#   define norm2 F77_FUNC(psnorm2, PSNORM2)
#   define REAL_FMT "f"
#else
    typedef double real;
#   define REAL_MPI_TYPE MPI_DOUBLE
#   define blas_dot F77_FUNC(ddot, DDOT)
#   define blas_gemv F77_FUNC(dgemv, DGEMV)
#   define blas_gemm F77_FUNC(dgemm, DGEMM)
#   define blas_axpy F77_FUNC(daxpy, DAXPY)
#   define lapack_syev F77_FUNC(dsyev, DSYEV)
#   define lapack_gesvd F77_FUNC(dgesvd, DGESVD)
#   define lapack_getrf F77_FUNC(dgetrf, DGETRF)
#   define lapack_getri F77_FUNC(dgetri, DGETRI)
#   define arpack_saupd F77_FUNC(dsaupd, DSAUPD)
#   define arpack_seupd F77_FUNC(dseupd, DSEUPD)
#   define parpack_psaupd F77_FUNC(pdsaupd, PDSAUPD)
#   define parpack_pseupd F77_FUNC(pdseupd, PDSEUPD)
#   define parpack_pnorm2 F77_FUNC(pdnorm2, PDNORM2)
#   define REAL_FMT "d"
#endif // OPSEC_SINGLE


/* C declarations of BLAS, LAPACK, and PARPACK Fortran routines. */

extern real blas_dot(int* n, real* x, int* incx, real* y, int* incy);

extern void blas_gemv(char* trans, int* m, int* n, real* alpha, real* a, int* lda,
                      real* x, int* incx, real* beta, real* y, int* incy);

extern void blas_gemm(char* transa, char* transb, int* m, int* n, int* k,
                      real* alpha, real* a, int* lda, real* b, int* ldb, real* beta, real* c, int* ldc);

extern void blas_axpy(int* nloc, real* a, real* x, int* incx, real* y, int* incy);

extern void lapack_syev(char* jobz, char* uplo, int* n, real* a, int* lda,
                        real* w, real* work, int* lwork, int* info);

extern void lapack_gesvd(char* jobu, char* jobvt, int* m, int* n, real* a, int* lda, real* s,
                         real* u, int* ldu, real* vt, int* lvt, real* work, int* lwork, int* info);

extern void lapack_getrf(int* m, int* n, real* a, int* lda, int* ipiv, int* info);

extern void lapack_getri(int* n, real* a, int* lda, int* ipiv, real* work, int* lwork, int* info);

extern void arpack_saupd(int *ido, char *bmat, int *n, char *which,
                         int *nev, real *tol, real *resid, int *ncv, real *V,
                         int *ldv, int *iparam, int *ipntr, real *workd,
                         real *workl, int *lworkl, int *info);

extern void arpack_seupd(int *rvec, char *howmny, int *select, real *d,
                         real *z, int *ldz, real *sigma, char *bmat, int *n,
                         char *which, int *nev, real *tol, real *resid, int *ncv,
                         real *v, int *ldv, int *iparam, int *ipntr, real *workd,
                         real *workl, int *lworkl, int *info);

extern void parpack_psaupd(int *comm, int *ido, char *bmat, int *n, char *which,
                           int *nev, real *tol, real *resid, int *ncv, real *V,
                           int *ldv, int *iparam, int *ipntr, real *workd,
                           real *workl, int *lworkl, int *info);

extern void parpack_pseupd(int *comm, int *rvec, char *howmny, int *select, real *d,
                           real *z, int *ldz, real *sigma, char *bmat, int *n,
                           char *which, int *nev, real *tol, real *resid, int *ncv,
                           real *v, int *ldv, int *iparam, int *ipntr, real *workd,
                           real *workl, int *lworkl, int *info);

extern real parpack_pnorm2(int* comm, int* nloc, real* x, int* incx);

#ifdef __cplusplus
}
#endif

#endif // OPSEC_H
