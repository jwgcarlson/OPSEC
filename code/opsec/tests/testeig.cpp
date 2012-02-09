/* testeig
 *
 * Test program for eigenvalue/eigenvector computation using PARPACK. */

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <mpi.h>

#include "CuTest.h"
#include "eig.h"
#include "opsec.h"

/* Test that matrix-vector multiplication using BLAS works properly */
void TestGemv(CuTest* tc) {
    char trans = 'T';
    int m = 4;
    int n = 3;
    real alpha = 1.0;
    real beta = 0.0;
    int incx = 1;
    int incy = 1;
    real A[3*4] = { 1, 2, 3, 4,
                    -1, 1, -1, 1,
                    0, 0, 0, 1 };
    real x[4] = { 1, 2, 0, -1 };
    real y[3];
    blas_gemv(&trans, &m, &n, &alpha, A, &m, x, &incx, &beta, y, &incy);

    /* A*x should equal { 1, 0, -1 } */
    CuAssertDblEquals(tc, 1, y[0], 1e-10);
    CuAssertDblEquals(tc, 0, y[1], 1e-10);
    CuAssertDblEquals(tc, -1, y[2], 1e-10);
}

/* Test that matrix-matrix multiplication using BLAS works properly */
void TestGemm(CuTest* tc) {
    char trans = 'T', notrans = 'N';
    real alpha = 1.0, beta = 0.0;

    int m = 2;
    int n = 3;
    int k = 6;
    real B[2*6] = { 0.5, 0.5, 0.5, 0.5, 0, 0,
                    -0.5, 0, -0.5, 0, 0.5, 0.5 };
    real C[3*6] = { 1, 2, 3, 4, 5, 6,
                    -1, 1, -1, 1, -1, 1,
                    0, 0, 0, 0, 0, 1 };
    real tmp[3*2];
    real A[2*2];

//    printf("B = %g %g %g %g %g %g %g %g %g %g %g %g\n", B[0], B[1], B[2], B[3], B[4], B[5], B[6], B[7], B[8], B[9], B[10], B[11]);
    blas_gemm(&trans, &notrans, &n, &m, &k, &alpha, C, &k, B, &k, &beta, tmp, &n);
//    printf("tmp = %g %g %g %g %g %g\n", tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5]);

    real* Bloc = B + 0;
    blas_gemm(&trans, &notrans, &m, &m, &n, &alpha, tmp, &n, Bloc, &k, &beta, A, &m);

    /* A = Bloc * C * B^T should equal { 2.5, 2.5, -2.5, -2 } */
//    printf("A = %g %g %g %g\n", A[0], A[1], A[2], A[3]);
    CuAssertDblEquals(tc, 2.5, A[0], 1e-10);
    CuAssertDblEquals(tc, 2.5, A[1], 1e-10);
    CuAssertDblEquals(tc, -2.5, A[2], 1e-10);
    CuAssertDblEquals(tc, -2, A[3], 1e-10);

    Bloc = B + 3;
    blas_gemm(&trans, &notrans, &m, &m, &n, &alpha, tmp, &n, Bloc, &k, &beta, A, &m);
    /* A = Bloc * C * B^T should equal { 2.5, 1.75, 0, 0.75 } */
//    printf("A = %g %g %g %g\n", A[0], A[1], A[2], A[3]);
    CuAssertDblEquals(tc, 2.5, A[0], 1e-10);
    CuAssertDblEquals(tc, 1.75, A[1], 1e-10);
    CuAssertDblEquals(tc, 0, A[2], 1e-10);
    CuAssertDblEquals(tc, 0.75, A[3], 1e-10);
}

/* J_x |m> = (1/2) * [ sqrt{j(j+1) - m(m+1)} |m+1> + sqrt{j(j+1) - m(m-1)} |m-1> ] */
void TestEig(CuTest* tc) {
    int n = 11, nev = 6, ncv = 8, nconv;
    real *A, *evals, *evecs, *v;
    int a, b;
    real m, j = (n - 1.)/2.;

    A = (real*) opsec_malloc(n*n*sizeof(real));
    for(a = 0; a < n; a++) {
        for(b = 0; b < n; b++) {
            m = j - b;
            if(a == b-1)
                A[b + a*n] = 0.5 * sqrt(j*(j+1) - m*(m+1));
            else if(a == b+1)
                A[b + a*n] = 0.5 * sqrt(j*(j+1) - m*(m-1));
            else
                A[b + a*n] = 0;
        }
    }

    evals = (real*) opsec_malloc(nev*sizeof(real));
    evecs = (real*) opsec_malloc(n*ncv*sizeof(real));
    nconv = eig(n, nev, ncv, A, evals, evecs);

    printf("Asked for %d eigenvalues, got %d\n", nev, nconv);

    /* Print eigenvalues and eigenvectors */
    for(a = 0; a < nconv; a++) {
        printf("lambda_%d = %5.3f, v_%d = {", a, evals[a], a);
        for(b = 0; b < n; b++)
            printf(" %5.3f", evecs[a*n + b]);
        printf(" }\n");
    }

    CuAssertTrue(tc, fabs(evals[0] - 0.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[1] - 1.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[2] - 2.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[3] - 3.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[4] - 4.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[5] - 5.0) < 1e-10);

    free(A);
    free(evals);
    free(evecs);
}

void TestPeig(CuTest* tc) {
    int n = 11, nev = 6, ncv = 8, nloc, nconv;
    real *A, *evals, *evecs, *v;
    int a, aloc, b;
    int nprocs, myid;
    int *locsizes, *locdisps;
    real m, j = (n - 1.)/2.;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    locsizes = (int*) malloc(nprocs*sizeof(int));
    locdisps = (int*) malloc(nprocs*sizeof(int));
    for(int p = 0; p < nprocs; p++) {
        locsizes[p] = (n/nprocs) + (p < (n % nprocs));
        locdisps[p] = (p == 0) ? 0 : locdisps[p-1] + locsizes[p-1];
    }
    assert(locdisps[nprocs-1] + locsizes[nprocs-1] == n);
    nloc = locsizes[myid];

    A = (real*) opsec_malloc(n*nloc*sizeof(real));
    for(aloc = 0; aloc < nloc; aloc++) {
        a = aloc + locdisps[myid];
        for(b = 0; b < n; b++) {
            m = j - b;
            if(a == b-1)
                A[b + aloc*n] = 0.5 * sqrt(j*(j+1) - m*(m+1));
            else if(a == b+1)
                A[b + aloc*n] = 0.5 * sqrt(j*(j+1) - m*(m-1));
            else
                A[b + aloc*n] = 0;
        }
    }

    evals = (real*) opsec_malloc(nev*sizeof(real));
    evecs = (real*) opsec_malloc(nloc*ncv*sizeof(real));
    nconv = peig(MPI_COMM_WORLD, n, nloc, nev, ncv, A, evals, evecs);

    if(myid == 0) {
        printf("Asked for %d eigenvalues, got %d\n", nev, nconv);
        v = (real*) opsec_malloc(n*sizeof(real));
    }

    /* Print eigenvalues and eigenvectors */
    for(a = 0; a < nconv; a++) {
        MPI_Gatherv(evecs + a*nloc, nloc, REAL_MPI_TYPE, v, locsizes, locdisps, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
        if(myid == 0) {
            printf("lambda_%d = %5.3f, v_%d = {", a, evals[a], a);
            for(b = 0; b < n; b++)
                printf(" %5.3f", v[b]);
            printf(" }\n");
        }
    }

    CuAssertTrue(tc, fabs(evals[0] - 0.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[1] - 1.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[2] - 2.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[3] - 3.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[4] - 4.0) < 1e-10);
    CuAssertTrue(tc, fabs(evals[5] - 5.0) < 1e-10);

    free(locsizes);
    free(locdisps);
    free(A);
    free(evals);
    free(evecs);
    if(myid == 0)
        free(v);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestGemv);
    SUITE_ADD_TEST(suite, TestEig);
    SUITE_ADD_TEST(suite, TestPeig);
    SUITE_ADD_TEST(suite, TestGemm);

    return suite;
}

int main(int argc, char **argv) {
    int myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    CuString* output = CuStringNew();
    CuSuite* suite = CuSuiteNew();

    CuSuiteAddSuite(suite, GetSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    if(myid == 0)
        printf("%s\n", output->buffer);

    MPI_Finalize();
    return (suite->failCount == 0) ? 0 : 1;
}
