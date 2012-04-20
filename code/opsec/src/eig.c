#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "abn.h"
#include "cfg.h"
#include "eig.h"

/* Serial eigenvalue routine.  Uses ARPACK. */
int eig(int n, int nev, int ncv, real *A, real *evals, real *evecs) {
    int iparam[11], ipntr[11], ierr;
    real *workl, *workd, *resid, *xfull;
    int ido, lworkl, info;
    real tol, sigma;
    real *x, *y;
    char bmat[] = "I";
    char which[] = "LA";
    int rvec = 1;
    char howmny[] = "All";
    int *select;
    int itr;

    /* Set PARPACK parameters */
    tol = 1e-5;
    iparam[0] = 1;      // ishfts
    iparam[2] = 300;    // maxitr
    iparam[6] = 1;      // mode
    lworkl = ncv*(ncv + 8);
    workl = (real*) opsec_calloc(lworkl, sizeof(real));
    workd = (real*) opsec_calloc(3*n, sizeof(real));
    resid = (real*) opsec_calloc(n, sizeof(real));
    xfull = (real*) opsec_calloc(n, sizeof(real));
    select = (int*) opsec_calloc(ncv, sizeof(int));

    /* Begin reverse communication loop */
    itr = 0;
    info = 0;
    ido = 0;
    while(ido != 99) {
        arpack_saupd(&ido, bmat, &n, which, &nev,
                     &tol, resid, &ncv, evecs, &n, iparam, ipntr,
                     workd, workl, &lworkl, &info);

        if(ido == 1 || ido == -1) {
            /* GEMV parameters.  Since Fortran expects matrices to be stored in
             * column-major order (while ours uses row-major), we treat A as an
             * n-by-nloc matrix and compute the transpose product A'*x. */
            static char trans = 'T';
            static int incx = 1;
            static int incy = 1;
            static real alpha = 1.0;
            static real beta = 0.0;

            /* Compute y = A*x (and don't forget Fortran indexing conventions!) */
            x = &workd[ipntr[0] - 1];
            y = &workd[ipntr[1] - 1];

//            printf("Iteration %d\n", ++itr);

            blas_gemv(&trans, &n, &n, &alpha, A, &n, x, &incx, &beta, y, &incy);
        }
    }

    printf("Number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);

    /* Check return code */
    if(info < 0) {
        /* Error encountered.  Abort. */
        fprintf(stderr, "arpack_saupd returned error: info = %d\n", info);
        return info;
    }
    else {
        arpack_seupd(&rvec, howmny, select, evals, evecs, &n, &sigma,
                     bmat, &n, which, &nev, &tol, resid, &ncv, evecs, &n,
                     iparam, ipntr, workd, workl, &lworkl, &ierr);

        if(ierr != 0)
            fprintf(stderr, "arpack_seupd returned error: ierr = %d\n", ierr);
    }

#if 0
    {
        /* Debugging: check residuals  || A*x - lambda*x || */
        y = (real*) opsec_calloc(n, sizeof(real));
        for(int i = 0; i < iparam[4]; i++) { 
            char trans = 'T';
            int incx = 1;
            int incy = 1;
            real alpha = 1.0;
            real beta = 0.0;
            real a = -evals[i];
            real* x = (real*) opsec_calloc(n, sizeof(real));
            gemv(&trans, &n, &n, &alpha, A, &n, x, &incx, &beta, y, &incy);
            axpy(&nloc, &a, evecs + i*nloc, &incx, y, &incy);
            real d = norm2(&fcomm, &nloc, y, &incy);
            if(myid == 0)
                printf("Eigenvalue %d: lambda = %16.16f, |A*x - lambda*x| = %16.16f\n", i, evals[i], d);
            ierr = MPI_Allgatherv(y, nloc, REAL_MPI_TYPE, xfull, locsizes, locdisps, REAL_MPI_TYPE, comm);
            if(myid == 0) {
                printf("y =");
                for(int i = 0; i < n; i++)
                    printf(" %g", xfull[i]);
                printf("\n");
            }
        }
        free(y);
    }
#endif

    free(workl);
    free(workd);
    free(resid);
    free(xfull);
    free(select);

    return iparam[4];
}

#ifdef OPSEC_USE_MPI

/* Parallel eigenvalue routine.  Uses MPI and PARPACK. */
int peig(int n, int nloc, int nev, int ncv, real tol, real *A, real *evals, real *evecs) {
    MPI_Fint fcomm;
    int myid, nprocs;
    int *locsizes, *locdisps;
    int iparam[11], ipntr[11], ierr;
    real *workl, *workd, *resid, *xfull;
    int ido, lworkl, info;
    real sigma;
    real *x, *y;
    char bmat[] = "I";
    char which[] = "LA";
    int rvec = 1;
    char howmny[] = "All";
    int *select;
    int p, itr;

    /* Set PARPACK parameters */
    iparam[0] = 1;      // ishfts
    iparam[2] = 1000;   // maxitr
    iparam[6] = 1;      // mode
    lworkl = ncv*(ncv + 8);
    workl = (real*) opsec_calloc(lworkl, sizeof(real));
    workd = (real*) opsec_calloc(3*nloc, sizeof(real));
    resid = (real*) opsec_calloc(nloc, sizeof(real));
    xfull = (real*) opsec_calloc(n, sizeof(real));
    select = (int*) opsec_calloc(ncv, sizeof(int));

    /* Get MPI info */
    fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Determine local problem lengths, by gathering nloc from each process */
    locsizes = (int*) opsec_malloc(nprocs*sizeof(int));
    locdisps = (int*) opsec_malloc(nprocs*sizeof(int));
    MPI_Allgather(&nloc, 1, MPI_INT, locsizes, 1, MPI_INT, MPI_COMM_WORLD);
    for(p = 0; p < nprocs; p++)
        locdisps[p] = (p == 0) ? 0 : locdisps[p-1] + locsizes[p-1];
    assert(locdisps[nprocs-1] + locsizes[nprocs-1] == n);

    /* Begin reverse communication loop */
    itr = 0;
    info = 0;
    ido = 0;
    while(ido != 99) {
        parpack_psaupd(&fcomm, &ido, bmat, &nloc, which, &nev,
                       &tol, resid, &ncv, evecs, &nloc, iparam, ipntr,
                       workd, workl, &lworkl, &info);

        if(ido == 1 || ido == -1) {
            /* GEMV parameters.  Since Fortran expects matrices to be stored in
             * column-major order (while ours uses row-major), we treat A as an
             * n-by-nloc matrix and compute the transpose product A'*x. */
            static char trans = 'T';
            static int incx = 1;
            static int incy = 1;
            static real alpha = 1.0;
            static real beta = 0.0;

            /* Compute y = A*x (and don't forget Fortran indexing conventions!) */
            x = &workd[ipntr[0] - 1];
            y = &workd[ipntr[1] - 1];

//            if(myid == 0) {
//                printf("pdsaupd call number %d\n", ++itr);
//                fflush(stdout);
//            }

            MPI_Allgatherv(x, nloc, REAL_MPI_TYPE, xfull, locsizes, locdisps, REAL_MPI_TYPE, MPI_COMM_WORLD);
            blas_gemv(&trans, &n, &nloc, &alpha, A, &n, xfull, &incx, &beta, y, &incy);
        }
    }

    if(myid == 0) {
        printf("Number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
        printf("  info = %d\n", info);
        printf("  nconv = %d, nev = %d\n", iparam[4], nev);
    }

    if(myid == 0) {
        time_t t = time(NULL);
        printf("Time: %s\n", ctime(&t));
        printf("Post-processing Ritz values and vectors\n");
    }

    /* Check return code */
    if(info < 0) {
        /* Error encountered.  Abort. */
        if(myid == 0)
            fprintf(stderr, "parpack_psaupd returned error: info = %d\n", info);
        return info;
    }
    else {
        parpack_pseupd(&fcomm, &rvec, howmny, select, evals, evecs, &nloc, &sigma,
                       bmat, &nloc, which, &nev, &tol, resid, &ncv, evecs, &nloc,
                       iparam, ipntr, workd, workl, &lworkl, &ierr);

        if(ierr != 0) {
            if(myid == 0)
                fprintf(stderr, "parpack_pseupd returned error: ierr = %d\n", ierr);
        }
    }

    if(myid == 0) {
        time_t t = time(NULL);
        printf("Time: %s\n", ctime(&t));
    }

#if 0
    {
        int i;
        /* Debugging: check residuals  || A*x - lambda*x || */
        y = (real*) opsec_calloc(nloc, sizeof(real));
        for(i = iparam[4]-1; i >= 0; i--) { 
            static char trans = 'T';
            static int incx = 1;
            static int incy = 1;
            static real alpha = 1.0;
            static real beta = 0.0;
            real a = -evals[i];
            ierr = MPI_Allgatherv(&evecs[i*nloc], nloc, REAL_MPI_TYPE, xfull, locsizes, locdisps, REAL_MPI_TYPE, MPI_COMM_WORLD);
            blas_gemv(&trans, &n, &nloc, &alpha, A, &n, xfull, &incx, &beta, y, &incy);
            blas_axpy(&nloc, &a, &evecs[i*nloc], &incx, y, &incy);
            real d = parpack_pnorm2(&fcomm, &nloc, y, &incy);
            if(myid == 0)
                printf("Eigenvalue %d: lambda = %16.16f, |A*x - lambda*x| = %16.16f\n", iparam[4]-i, evals[i], d);
            ierr = MPI_Allgatherv(y, nloc, REAL_MPI_TYPE, xfull, locsizes, locdisps, REAL_MPI_TYPE, MPI_COMM_WORLD);
        }
        free(y);
    }
#endif

    /* Clean up */
    free(locsizes);
    free(locdisps);
    free(workl);
    free(workd);
    free(resid);
    free(xfull);
    free(select);

    return iparam[4];
}

#endif // OPSEC_USE_MPI
