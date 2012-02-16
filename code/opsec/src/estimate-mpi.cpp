/* estimate-mpi
 *
 * Produce parameter estimates $\hat{p}_m$ from pixel values, as well as
 * a covariance matrix and various other quantities.
 *
 * This program is designed to be run in parallel using MPI. */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "Matrix.h"
#include "Model.h"
#include "SplitFile.h"

const char* usage =
    "Usage: estimate [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h             Display this help\n"
    "  -c FILE        Read additional configuration options from FILE\n"
    "Configuration options:\n"
    "  estfile=FILE   Write parameter estimates to FILE (required)\n"
    "  evalfile=FILE  Read S/N eigenvalues from FILE (required)\n"
    "  covfile=FILE   Read covariance matrix derivatives from FILE (required)\n"
    "  pixfile=FILE   Read pixel values from FILE (required)\n"
    "  mixing=MODE    Mixing matrix choice: identity, inverse (default), or inverse sqrt\n";

/* Read signal-to-noise eigenvalues from file on root process, then broadcast
 * to all processes. */
real* ReadEvals(const char* evalfile, int* Nmodes_ = NULL) {
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    int Nmodes = -1;
    real* evals = NULL;
    if(me == 0) {
        FILE* feval = fopen(evalfile, "rb");
        if(feval == NULL) {
            fprintf(stderr, "ReadEvals: could not read from file '%s'\n", evalfile);
        }
        else {
            /* File successfully opened */
            Config opts = cfg_new();
            size_t n, size;
            int err = abn_read(feval, (void**)&evals, &n, &size, NULL, NULL, opts);
            if(err || evals == NULL) {
                fprintf(stderr, "ReadEvals: error reading signal-to-noise eigenvalues from '%s'\n", evalfile);
            }
            else {
                /* Eigenvalues successfully read */
                int Nmodes = cfg_get_int(opts, "Nmodes");
                if((int)n != Nmodes)
                    fprintf(stderr, "ReadEvals: sanity check fail: n = %zd, Nmodes = %d\n", n, Nmodes);
                if(size != sizeof(real))
                    fprintf(stderr, "ReadEvals: sanity check fail: size = %zd, sizeof(real) = %zd\n", size, sizeof(real));
            }
            fclose(feval);
            cfg_destroy(opts);
        }
    }

#ifdef OPSEC_USE_MPI
    /* Broadcast Nmodes and evals to all processes */
    MPI_Bcast(&Nmodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(Nmodes >= 0) {
        if(me != 0)
            evals = (real*) opsec_malloc(Nmodes*sizeof(real));
        MPI_Bcast(evals, Nmodes, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
    }
#endif

    if(Nmodes_ != NULL)
        *Nmodes_ = Nmodes;
    return evals;
}

/* Read pixel values from file. */
real* ReadPixels(const char* pixfile, int* Npixels_ = NULL) {
    FILE* fpix = fopen(pixfile, "rb");
    if(fpix == NULL) {
        fprintf(stderr, "ReadPixels: could not read from file '%s'\n", pixfile);
        return NULL;
    }

    size_t n, size;
    real* pixels = NULL;
    int err = abn_read(fpix, (void**)&pixels, &n, &size, NULL, NULL, NULL);
    if(err || pixels == NULL) {
        fprintf(stderr, "ReadPixels: error reading from '%s'\n", pixfile);
    }
    else {
        if(size != sizeof(real))
            fprintf(stderr, "ReadPixels: sanity check fail: size = %zd, sizeof(real) = %zd\n", size, sizeof(real));

        if(Npixels_)
            *Npixels_ = (int)n;
    }

    fclose(fpix);
    return pixels;
}

/* Read an m-by-n matrix from file, and distribute across processes.  The
 * matrix is distributed by rows.  The local matrix A must have enough memory
 * to store the local portion of the matrix. */
template<class ScalarType>
void ReadMatrix(SplitFile& f, int m, int n, MyMatrix<ScalarType>& A, int nprocs = 1, int me = 0) {
#ifdef OPSEC_DEBUG
    assert(A.n == n && A.m == (m/nprocs) + (me < (m % nprocs)));
#endif

    /* Compute number of rows owned by each process, and their offsets */
    std::vector<int> mloc(nprocs), offset(nprocs);
    for(int p = 0; p < nprocs; p++) {
        mloc[p] = (m/nprocs) + (p < (m % nprocs));
        offset[p] = p*(m/nprocs) + std::min(p, m % nprocs);
    }

    if(me == 0) {
        for(int p = 0; p < nprocs; p++) {
        }
    }
}

/* Read in derivatives of the covariance matrix from file. */
void ReadDerivatives(const char* covfile, bool multifile, int* Nparams_ = NULL, int* Nmodes_ = NULL) {
    int Nparams = 0;
    int Nmodes = 0;
    real** Cp = NULL;
    FILE* fcov = NULL;
    Config opts = NULL;
    if(!multifile) {
        /* Single file */

        /* Open file */
        fcov = fopen(covfile, "rb");
        if(fcov == NULL) {
            fprintf(stderr, "ReadDerivatives: could not read from file '%s'\n", covfile);
            return NULL;
        }

        /* Read header */
        opts = cfg_new();
        size_t n, size;
        if(abn_read_header(fcov, &n, &size, NULL, NULL, opts) != 0) {
            fprintf(stderr, "ReadDerivatives: error reading from '%s'\n", covfile);
            return NULL;
        }

        /* Set parameters, allocate memory, and read in data */
        Nparams = cfg_get_int(opts, "Nparams");
        Nmodes = cfg_get_int(opts, "Nmodes");
        Cp = (real**) opsec_malloc(Nparams*sizeof(real*));
        for(int m = 0; m < Nparams; m++) {
            Cp[m] = (real*) opsec_malloc(Nmodes*Nmodes*sizeof(real));
            size_t nread = fread(Cp[m], sizeof(real), Nmodes*Nmodes, fcov);
            if(nread != Nmodes*Nmodes) {
                fprintf(stderr, "ReadDerivatives: error reading Cp[%d] from '%s'\n", m, covfile);
                return NULL;
            }
        }

        fclose(fcov);
        cfg_destroy(opts);
    }
    else {
        /* Multi-file */

        char filename[512];

        /* Loop through each file */
        Nparams = 1;    // this should be set properly after reading the first file
        for(int m = 0; m < Nparams; m++) {
            /* Open file */
            snprintf(filename, sizeof(filename), "%s.%03d", covfile, m);
            fcov = fopen(filename, "rb");
            if(fcov == NULL) {
                fprintf(stderr, "ReadDerivatives: could not read from '%s'\n", filename);
                return NULL;
            }

            /* Read header */
            opts = cfg_new();
            size_t n, size;
            if(abn_read_header(fcov, &n, &size, NULL, NULL, opts) != 0) {
                fprintf(stderr, "ReadDerivatives: error reading from '%s'\n", filename);
                return NULL;
            }

            if(m == 0) {
                /* Use the header of the first file to set parameters */
                Nparams = cfg_get_int(opts, "Nparams");
                Nmodes = cfg_get_int(opts, "Nmodes");
                assert(Nparams > 0 && Nmodes > 0);

                /* TODO: check memory availibility, to make sure there will be no swapping */

                /* Allocate memory for the matrices */
                Cp = (real**) opsec_malloc(Nparams*sizeof(real*));
                for(int m = 0; m < Nparams; m++)
                    Cp[m] = (real*) opsec_malloc(Nmodes*Nmodes*sizeof(real));
            }
            else {
                /* Make sure the parameters for later files match the first file */
                if(cfg_has_key(opts, "Nparams"))
                    assert(cfg_get_int(opts, "Nparams") == Nparams);
                if(cfg_has_key(opts, "Nmodes"))
                    assert(cfg_get_int(opts, "Nmodes") == Nmodes);
            }

            /* Read in the data */
            size_t nread = fread(Cp[m], sizeof(real), Nmodes*Nmodes, fcov);
            if(nread != Nmodes*Nmodes) {
                fprintf(stderr, "ReadDerivatives: error reading Cp[%d] from '%s'\n", m, filename);
                return NULL;
            }

            cfg_destroy(opts);
            fclose(fcov);
        }

    }

    if(Nparams_)
        *Nparams_ = Nparams;
    if(Nmodes_)
        *Nmodes_ = Nmodes;

    return Cp;
}

void matrix_print(const Matrix<real>& A) {
    for(int i = 0; i < A.M; i++) {
        printf("%s[ ", (i == 0) ? "[" : " ");
        for(int j = 0; j < A.N; j++)
            printf("% 9.5e ", A(i,j));
        printf("]%s\n", (i == A.M-1) ? "]" : "");
    }
}

/* Return the square root of the symmetric, positive-definite, N-by-N matrix A. */
Matrix<real> matrix_sqrt(const Matrix<real>& A) {
    int N = A.M;
    assert(A.N == N);
    /* Eigenvalues (output only). */
    vector<real> w(N);
    /* On input, the symmetric matrix A.  On output, the orthonormal
     * eigenvectors of A, stored in column-major format. */
    Matrix<real> Ut = A;
    int lwork = 4*N;
    real* work = (real*)opsec_malloc(lwork*sizeof(real));
    int info;
    char jobz = 'V', uplo = 'U';
    lapack_syev(&jobz, &uplo, &N, &Ut(0,0), &N, &w[0], work, &lwork, &info);
    if(info != 0) {
        fprintf(stderr, "matrix_sqrt: lapack_syev returned %d\n", info);
    }
    Matrix<real> U = transpose(Ut);
    
    /* Multiply U by sqrt(S). */
    int i, j;
    #pragma omp parallel for
    for(i = 0; i < N; i++)
        w[i] = sqrt(w[i]);
    #pragma omp parallel for private(j)
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            U(i,j) *= w[j];

    /* Multiply (U*sqrt(S)) by U^T. */
    return dot(U, Ut);
}

/* Return the inverse of the N-by-N matrix A. */
Matrix<real> matrix_inv(const Matrix<real>& A) {
    int N = A.M;
    assert(A.N == N);

    /* On input, the matrix A.  On output, the factors L and U from the LU
     * factorization A = P*L*U. */
    real* lu = (real*)opsec_malloc(N*N*sizeof(real));
    memcpy(lu, &A(0,0), N*N*sizeof(real));

    /* On output, the pivot indices that determine the permutation matrix P. */
    int* ipiv = (int*)opsec_malloc(N*sizeof(int));

    /* Perform LU factorization on the matrix A. */
    int info;
    lapack_getrf(&N, &N, lu, &N, ipiv, &info);
    if(info != 0) {
        fprintf(stderr, "matrix_inv: lapack_getrf returned %d\n", info);
    }

    /* Use the LU factorization to compute the inverse of A. */
    int lwork = 4*N;
    real* work = (real*)opsec_malloc(lwork*sizeof(real));
    lapack_getri(&N, lu, &N, ipiv, work, &lwork, &info);
    if(info != 0) {
        fprintf(stderr, "matrix_inv: lapack_getri returned %d\n", info);
    }

    free(ipiv);
    free(work);
    return Matrix<real>(N, N, lu, false);
}

int main(int argc, char* argv[]) {
    /* Initialize MPI */
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* Read configuration options */
    Config cfg = opsec_init(argc, argv, usage);

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "estfile,evalfile,covfile,pixfile,Nmodes", ",")) {
        fprintf(stderr, "estimate: missing configuration options\n");
        fputs(usage, stderr);
        opsec_exit(1);
    }

    const char* estfile = cfg_get(cfg, "estfile");
    const char* evalfile = cfg_get(cfg, "evalfile");
    const char* covfile = cfg_get(cfg, "covfile");
    const char* pixfile = cfg_get(cfg, "pixfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");

    int m, n, i, j;

    string mixing = "inverse";
    if(cfg_has_key(cfg, "mixing"))
        mixing = cfg_get(cfg, "mixing");

    /* Determine output format for derivatives; assume single file by default */
    string comma_output = cfg_has_key(cfg, "comma_output") ? cfg_get(cfg, "comma_output") : "single";
    bool multifile = (comma_output.substr(0,5) == "multi");

    /* Load model, for prior parameter values */
    Model* model = InitializeModel(cfg);
    int Nparams = model->NumParams();


    int m = Nmodes;
    int n = Nparams;

    /* Determine the problem size and offset for each process */
    std::vector<int> mloc(nprocs), offset(nprocs);
    for(int p = 0; p < nprocs; p++) {
        nloc[p] = (m/nprocs) + (p < (m % nprocs));
        offset[p] = p*(m/nprocs) + std::min(p, m % nprocs);
    }

    /* Allocate local storage for matrices and vectors */
    real* svalues = NULL;       // signal-to-noise eigenvalues, allocated by ReadEvals()
    real* dvalues = (real*) opsec_malloc(n*m*mloc[me]*sizeof(real));    // covariance matrix derivatives
    real* rvalues = (real*) opsec_malloc(m*mloc[me]*sizeof(real));      // residue matrix

    /* Read signal-to-noise eigenvalues $s_i$ */
    int Nmodes_check;
    svalues = ReadEvals(evalfile, &Nmodes_check);
    if(svalues == NULL || Nmodes_check != Nmodes) {
        fprintf(stderr, "estimate-mpi: error reading signal-to-noise eigenvalues\n");
        opsec_exit(1);
    }

    /* Read covariance matrix derivatives $C_{,m}$ */
    int Nparams_check;
    real** D = ReadDerivatives(covfile, multifile, &Nparams_check, &Nmodes_check);
    if(D == NULL || Nparams_check != Nparams || Nmodes_check != Nmodes) {
        fprintf(stderr, "estimate: error reading covariance derivatives\n");
        opsec_exit(1);
    }

    /* Compute diagonal elements of covariance matrix $C_{ii} = 1 + \lambda_i$ */
    for(i = 0; i < Nmodes; i++)
        c[i] = 1. + c[i];

    /* Read pixel values $y_i$ */
    int Npixels;
    real* y = ReadPixels(pixfile, &Npixels);
    if(y == NULL || Npixels != Nmodes) {
        fprintf(stderr, "estimate: error reading pixel values\n");
        opsec_exit(1);
    }

    real pn, *Dm, *Dn;

    /* Compute residue matrix $R_{ij} = C_{ij} - \sum_n C_{ij,n} p_n$ */
    real* R = (real*) calloc(Nmodes*Nmodes, sizeof(real));      // initialize to 0's
    for(i = 0; i < Nmodes; i++)
        R[i*Nmodes + i] = c[i];
    for(n = 0; n < Nparams; n++) {
        Dn = D[n];
        pn = model->GetParam(n);
        #pragma omp parallel for private(j)
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                R[i*Nmodes + j] -= Dn[i*Nmodes + j] * pn;
    }

    /* Pre-multiply covariance matrix derivatives by $C^{-1}$ */
    #pragma omp parallel for private(i,j,Dn)
    for(n = 0; n < Nparams; n++) {
        Dn = D[n];
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                Dn[i*Nmodes + j] /= c[i];
    }

    /* Compute Fisher matrix $F_{mn} = (1/2) \Tr[C^{-1} C_{,m} C^{-1} C_{,n}]$ */
    Matrix<real> F(Nparams, Nparams);
    #pragma omp parallel for private(m,n,i,j,Dm,Dn)
    for(int mn = 0; mn < Nparams*Nparams; mn++) {
        m = mn / Nparams;
        n = mn % Nparams;
        if(n > m)
            continue;
        Dm = D[m];
        Dn = D[n];
        F(m,n) = 0;
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                F(m,n) += Dm[i*Nmodes + j] * Dn[j*Nmodes + i];
        F(m,n) *= 0.5;
    }
    for(m = 0; m < Nparams; m++)
        for(n = m+1; n < Nparams; n++)
            F(m,n) = F(n,m);

    /* Post-multiply $D_n = C^{-1} C_{,n}$ by $C^{-1}$ */
    #pragma omp parallel for private(i,j,Dn)
    for(n = 0; n < Nparams; n++) {
        Dn = D[n];
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                Dn[i*Nmodes + j] /= c[j];
    }

    /* Compute quadratic estimators $q_n = (1/2) y^T C^{-1} C_{,n} C^{-1} y$ */
    Matrix<real> q(Nparams, 1);
    #pragma omp parallel for private(i,j,Dn)
    for(n = 0; n < Nparams; n++) {
        q(n) = 0;
        Dn = D[n];
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                q(n) += y[i] * Dn[i*Nmodes + j] * y[j];
        q(n) *= 0.5;
    }

    /* Compute bias corrections $f_n = (1/2) Tr[C^{-1} C_{,n} C^{-1} R]$ */
    Matrix<real> f(Nparams, 1);
    #pragma omp parallel for private(i,j,Dn)
    for(n = 0; n < Nparams; n++) {
        f(n) = 0;
        Dn = D[n];
        for(i = 0; i < Nmodes; i++)
            for(j = 0; j < Nmodes; j++)
                f(n) += Dn[i*Nmodes + j] * R[j*Nmodes + i];
        f(n) *= 0.5;
    }

    /* Free up memory */
    free(c);
    c = NULL;
    for(m = 0; m < Nparams; m++)
        free(D[m]);
    free(D);
    D = NULL;

    /* Choose mixing matrix */
    Matrix<real> M;
    if(mixing == "inverse sqrt")
        M = matrix_inv(matrix_sqrt(F));
    else if(mixing == "inverse")
        M = matrix_inv(F);
    else if(mixing == "identity")
        M = identity<real>(Nparams);
    else {
        fprintf(stderr, "estimate: unrecognized option: mixing = '%s'\n", mixing.c_str());
        fprintf(stderr, "estimate: using inverse\n");
        M = matrix_inv(F);
    }

    /* Compute window matrix */
    Matrix<real> W = dot(M, F);

    /* ... and normalize it */
    real rowsum;
    #pragma omp parallel for private(n,rowsum)
    for(m = 0; m < Nparams; m++) {
        rowsum = 0;
        for(n = 0; n < Nparams; n++)
            rowsum += W(m,n);
        for(n = 0; n < Nparams; n++)
            M(m,n) /= rowsum;
    }
    W = dot(M, F);

    /* Compute power estimates */
    Matrix<real> phat = dot(M, q - f);
    Matrix<real> cov = dot(M, dot(F, transpose(M)));

    /* Print them out */
    printf("p_n          error\n");
    for(n = 0; n < Nparams; n++)
        printf("%e %e\n", phat(n), sqrt(cov(n,n)));

    /* Write everything to binary file */
    FILE* fest = fopen(estfile, "wb");
    if(fest == NULL) {
        fprintf(stderr, "estimate: could not open '%s' for writing\n", estfile);
        opsec_exit(1);
    }

    fprintf(fest, "# Parameter estimates and associated matrices.\n");
    fprintf(fest, "# This file consists of several consecutive real-valued binary arrays:\n");
    fprintf(fest, "# - parameter estimates $\\hat{p}_n$ (vector of length Nparams)\n");
    fprintf(fest, "# - covariance matrix $Cov[\\hat{p}_m,\\hat{p}_n]$ (Nparams-by-Nparams matrix)\n");
    fprintf(fest, "# - Fisher matrix $F_{mn}$ (Nparams-by-Nparams matrix)\n");
    fprintf(fest, "# - window matrix $W_{mn}$ (Nparams-by-Nparams matrix)\n");
    fprintf(fest, "# - mixing matrix $M_{mn}$ (Nparams-by-Nparams matrix)\n");
    fprintf(fest, "# - quadratic combinations $q_n$ (vector of length Nparams)\n");
    fprintf(fest, "# - bias corrections $f_n$ (vector of length Nparams)\n");

    Config opts = cfg_new();
    cfg_set_int(opts, "Nparams", Nparams);

    abn_write(fest, &phat[0], Nparams, REAL_FMT, opts);
    abn_write(fest, &cov[0], Nparams*Nparams, REAL_FMT, NULL);
    abn_write(fest, &F[0], Nparams*Nparams, REAL_FMT, NULL);
    abn_write(fest, &W[0], Nparams*Nparams, REAL_FMT, NULL);
    abn_write(fest, &M[0], Nparams*Nparams, REAL_FMT, NULL);
    abn_write(fest, &q[0], Nparams, REAL_FMT, NULL);
    abn_write(fest, &f[0], Nparams, REAL_FMT, NULL);

    cfg_destroy(opts);
    fclose(fest);

    return 0;
}
