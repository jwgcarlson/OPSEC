#include <algorithm>
#include <cassert>
#include <cstring>
#include <ctime>

#include <mpi.h>

#include "ParpackSolver.h"

/* Return n divided by p, rounded up to nearest integer */
static inline int divup(int n, int p) {
    return (n + p - 1)/p;
}

ParpackSolver::ParpackSolver(int n_, MatrixFactory& matfact) {
    n = n_;

    int nprocs, me;
    Cblacs_pinfo(&me, &nprocs);
    nloc = n/nprocs + (me < (n % nprocs));

    /* Initialize BLACS process grid */
    pcontext = new slp::Context(nprocs, 1);

    /* Initialize matrix descriptor */
    Adesc = pcontext->new_descriptor(n, n, divup(n,nprocs), n);
    assert(nloc == Adesc->num_local_rows());
    assert(n == Adesc->num_local_cols());

    /* Allocate local memory for matrix $A$ */
    Avalues = (real*) opsec_malloc(Adesc->local_size() * sizeof(real));

    /* Fill in local matrix values */
    int nrows = Adesc->num_local_rows();
    int ncols = Adesc->num_local_cols();
    std::vector<int> rows(nrows), cols(ncols);
    for(int i = 0; i < nrows; i++)
        rows[i] = Adesc->row_l2g(i);
    for(int j = 0; j < ncols; j++)
        cols[j] = Adesc->col_l2g(j);
    matfact.ComputeMatrixValues(rows, cols, Avalues, Adesc->lld);

    /* Set default PARPACK parameters */
    nconv = -1;
    tol = 1e-5;
    maxitr = 1000;
    ncv = -1;

    /* These are initialized upon calling Solve() */
    xdesc = NULL;
    Bdesc = NULL;
    Bvalues = NULL;
}

ParpackSolver::~ParpackSolver() {
    delete xdesc;

    delete Bdesc;
    free(Bvalues);

    delete Adesc;
    free(Avalues);

    pcontext->exit();
    delete pcontext;
}

int ParpackSolver::Solve(int nev) {
    /* Get MPI info */
    int nprocs, me;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

    /* Select number of working Ritz vectors */
    if(ncv == -1)
        ncv = 2*nev;
    ncv = std::min(ncv, n-1);

    /* Initialize matrix descriptors */
    xdesc = pcontext->new_descriptor(n, 1, divup(n,nprocs), 1);
    Bdesc = pcontext->new_descriptor(n, nev, divup(n,nprocs), nev);
    assert(nloc == Bdesc->num_local_rows() && nloc == xdesc->num_local_rows());
    assert(ncv == Bdesc->num_local_cols() && 1 == xdesc->num_local_cols());

    /* Allocate local memory for eigenvector matrix $B$ */
    Bvalues = (real*) opsec_malloc(Bdesc->local_size() * sizeof(real));

    real sigma;
    int iparam[11], ipntr[11];

    /* Set PARPACK parameters */
    char bmat[] = "I";
    char which[] = "LA";
    char howmny[] = "All";
    iparam[0] = 1;      // ishfts
    iparam[2] = maxitr; // maxitr
    iparam[6] = 1;      // mode

    /* Allocate working memory */
    int lworkl = ncv*(ncv + 8);
    real* workl = (real*) opsec_calloc(lworkl, sizeof(real));
    real* workd = (real*) opsec_calloc(3*nloc, sizeof(real));
    real* resid = (real*) opsec_calloc(nloc, sizeof(real));
    int* select = (int*) opsec_calloc(ncv, sizeof(int));

    /* Begin reverse communication loop */
    int itr = 0;
    int info = 0;
    int ido = 0;
    while(ido != 99) {
        parpack_psaupd(&fcomm, &ido, bmat, &nloc, which, &nev,
                       &tol, resid, &ncv, Bvalues, &nloc, iparam, ipntr,
                       workd, workl, &lworkl, &info);

        if(ido == 1 || ido == -1) {
            /* Compute y = A*x (don't forget Fortran indexing conventions!) */
            slp::Matrix<real> A(Adesc, Avalues);
            slp::Matrix<real> x(xdesc, &workd[ipntr[0] - 1]);
            slp::Matrix<real> y(xdesc, &workd[ipntr[1] - 1]);
            slp::multiply(A, x, y);
        }
    }

    if(me == 0) {
        opsec_info("Number of Implicit Arnoldi update iterations taken is %d\n", iparam[2]);
        opsec_info("  info = %d\n", info);
        opsec_info("  nconv = %d, nev = %d\n", iparam[4], nev);

        time_t t = time(NULL);
        opsec_info("Time: %s\n", ctime(&t));
        opsec_info("Post-processing Ritz values and vectors\n");
    }

    /* Check return code */
    if(info < 0) {
        /* Error encountered.  Abort. */
        if(me == 0)
            opsec_error("parpack_psaupd returned error: info = %d\n", info);
        return info;
    }
    else {
        /* Save number of successfully computed eigenvalues */
        nconv = iparam[4];
        evals.resize(nconv);

        /* Retrieve eigenvalues and eigenvectors */
        int rvec = 1;
        int ierr;
        parpack_pseupd(&fcomm, &rvec, howmny, select, &evals[0], Bvalues, &nloc, &sigma,
                       bmat, &nloc, which, &nev, &tol, resid, &ncv, Bvalues, &nloc,
                       iparam, ipntr, workd, workl, &lworkl, &ierr);

        if(ierr != 0) {
            if(me == 0)
                opsec_error("parpack_pseupd returned error: ierr = %d\n", ierr);
        }
    }

    if(me == 0) {
        time_t t = time(NULL);
        opsec_info("Time: %s\n", ctime(&t));
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

#if 0
    /* Sort from largest to smallest eigenvalue */
    for(int j = 0; j < nconv/2; j++) {
        std::swap(evals[j], evals[nconv-j-1]);
        memcpy(workd, &B(0,j), nloc*sizeof(real));
        memcpy(&B(0,j), &B(0,nconv-j-1), nloc*sizeof(real));
        memcpy(&B(0,nconv-j-1), workd, nloc*sizeof(real));
    }
#endif

    /* Clean up */
    free(workl);
    free(workd);
    free(resid);
    free(select);

    return nconv;
}

const slp::Context* ParpackSolver::GetContext() const {
    return pcontext;
}

slp::Matrix<real> ParpackSolver::GetMatrix() const {
    return slp::Matrix<real>(Adesc, Avalues);
}

/* PARPACK's eigenvalues and eigenvectors are stored in reverse order, i.e.
 * sorted from smallest to largest.  We want to return them in
 * largest-to-smallest order, hence the indexing below. */
real ParpackSolver::GetEigenvalue(int i) const {
    assert(nconv > 0 && 0 <= i && i < nconv);
    return evals[nconv-i-1];
}

slp::Matrix<real> ParpackSolver::GetEigenvector(int i) const {
    assert(nconv > 0 && 0 <= i && i < nconv);
    return slp::Matrix<real>(xdesc, &Bvalues[(nconv-i-1)*Bdesc->lld]);
}

#if 0
const std::vector<real>& ParpackSolver::GetEigenvalues() const {
    return evals;
}

slp::Matrix<real> ParpackSolver::GetEigenvectors() const {
    return B;
}
#endif
