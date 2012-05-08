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
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "Model.h"
#include "MyMatrix.hpp"
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
 * to all processes.  The array evals must be able to hold Nevals elements. */
void ReadEvals(SplitFile& f, int Nevals, real* evals) {
    int me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    const char* evalfile = f.get_filename().c_str();
    if(me == 0) {
        /* Open evalfile */
        FILE* feval = fopen(evalfile, "rb");
        if(feval == NULL) {
            fprintf(stderr, "ReadEvals: error opening evalfile (%s,%d)\n", evalfile, errno);
            opsec_abort(1);
        }

        /* Read header */
        size_t n, size;
        char endian, fmt[ABN_MAX_FORMAT_LENGTH];
        Config opts = cfg_new();
        int err = abn_read_header(feval, &n, &size, &endian, fmt, opts);
        if(err) {
            fprintf(stderr, "ReadEvals: error reading header (%s,%d)\n", evalfile, err);
            opsec_abort(1);
        }
        if(size != sizeof(real) || strcmp(fmt, REAL_FMT) != 0) {
            fprintf(stderr, "ReadEvals: inconsistent data types (%s,%zd,%zd,%s,%s)\n", evalfile, size, sizeof(real), fmt, REAL_FMT);
            opsec_abort(1);
        }
        if((int) n != Nevals) {
            fprintf(stderr, "ReadEvals: inconsistent evals length (%s,%zd,%d)\n", evalfile, n, Nevals);
            opsec_abort(1);
        }
        cfg_destroy(opts);

        /* Read eigenvalues */
        size_t nread = fread(evals, sizeof(real), Nevals, feval);
        if((int) nread != Nevals) {
            fprintf(stderr, "ReadEvals: error reading data (%s,%zd,%d)\n", evalfile, nread, Nevals);
            opsec_abort(1);
        }
        fclose(feval);
    }

    /* Broadcast to other processes */
#ifdef OPSEC_USE_MPI
    MPI_Bcast(evals, Nevals, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
#endif
}

/* Read pixel values from file, and broadcast to all processes. */
real* ReadPixels(const char* pixfile, int* Npixels_ = NULL) {
    int me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    int Npixels = -1;
    real* pixels = NULL;
    if(me == 0) {
        FILE* fpix = fopen(pixfile, "rb");
        if(fpix == NULL) {
            fprintf(stderr, "ReadPixels: could not read from file '%s'\n", pixfile);
        }
        else {
            /* File opened successfully */
            size_t n, size;
            int err = abn_read(fpix, (void**)&pixels, &n, &size, NULL, NULL, NULL);
            if(err || pixels == NULL) {
                fprintf(stderr, "ReadPixels: error reading pixels from '%s'\n", pixfile);
            }
            else {
                /* Pixels read successfully */
                Npixels = (int) n;
                if(size != sizeof(real))
                    fprintf(stderr, "ReadPixels: sanity check fail: size = %zd, sizeof(real) = %zd\n", size, sizeof(real));

            }
        }
        fclose(fpix);
    }

#ifdef OPSEC_USE_MPI
    /* Broadcast Npixels and pixels to all processes */
    MPI_Bcast(&Npixels, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(Npixels >= 0) {
        if(me != 0)
            pixels = (real*) opsec_malloc(Npixels*sizeof(real));
        MPI_Bcast(pixels, Npixels, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
    }
#endif

    if(Npixels_)
        *Npixels_ = Npixels;
    return pixels;
}

/* Read ABN header of covariance matrix derivatives file.  If specified, check
 * that Nparams and Nmodes have the expected values.  This function should be
 * run by the root process only. */
void ReadCovHeader(SplitFile& fcov, int Nparams = -1, int Nmodes = -1, int packed = -1) {
    const char* covfile = fcov.get_filename().c_str();

    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    Config opts = cfg_new();
    int err = abn_read_header(fcov, &n, &size, &endian, fmt, opts);
    if(err) {
        fprintf(stderr, "ReadCovHeader: error reading covfile header (%s,%d)\n", covfile, err);
        opsec_exit(1);
    }
    if(endian != abn_endian()) {
        fprintf(stderr, "ReadCovHeader: endianness mismatch (%s,%c,%c)\n", covfile, abn_endian(), endian);
        opsec_exit(1);
    }
    if(size != sizeof(real) || strcmp(fmt, REAL_FMT) != 0) {
        fprintf(stderr, "ReadCovHeader: data type inconsistent (%s,%zd,%zd,%s,%s)\n", covfile, sizeof(real), size, REAL_FMT, fmt);
        opsec_exit(1);
    }
    if(Nparams > 0 && cfg_has_key(opts, "Nparams") && cfg_get_int(opts, "Nparams") != Nparams) {
        fprintf(stderr, "ReadCovHeader: Nparams inconsistent (%s,%d,%d)\n", covfile, Nparams, cfg_get_int(opts, "Nparams"));
        opsec_exit(1);
    }
    if(Nmodes > 0 && cfg_has_key(opts, "Nmodes") && cfg_get_int(opts, "Nmodes") != Nmodes) {
        fprintf(stderr, "ReadCovHeader: Nmodes inconsistent (%s,%d,%d)\n", covfile, Nmodes, cfg_get_int(opts, "Nmodes"));
        opsec_exit(1);
    }
    if(packed > 0 && cfg_has_key(opts, "packed") && cfg_get_int(opts, "packed") != packed) {
        fprintf(stderr, "ReadCovHeader: packed inconsistent (%s,%d,%d)\n", covfile, packed, cfg_get_int(opts, "packed"));
        opsec_exit(1);
    }
    cfg_destroy(opts);
}

/* Read a symmetric n-by-n matrix from file, and distribute across processes.
 * The matrix may either be stored in its entirety (all n*n elements), or in
 * column-major upper packed format, where only the n*(n+1)/2 unique elements
 * are stored.  Regardless, only the n*(n+1)/2 unique elements are returned,
 * distributed among processes as a naturally partitioned flat array. */
void ReadSymmetricMatrix(SplitFile& f, int n, real* ap, int packed) {
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* Total number of unique elements in a symmetric matrix */
    int ktotal = (n*(n+1))/2;

    /* Number of elements that will be assigned to each process (natural partition) */
    std::vector<int> kloc(nprocs);
    for(int p = 0; p < nprocs; p++)
        kloc[p] = (ktotal/nprocs) + (p < (ktotal % nprocs));

    /* Read in matrix elements, and distribute to the appropriate processes */
    if(me == 0) {
        /* Temporary storage for other process' elements.  We need room for at
         * most kloc[0] elements, since kloc[0] >= kloc[p] for all p >= 1.  We
         * only need this memory if nprocs > 1. */
        real* tmp = (nprocs == 1) ? NULL
                                  : (real*) opsec_malloc(kloc[0]*sizeof(real));

        /* Current column and row number, for non-packed storage */
        int col = 0, row = 0;

        for(int p = 0; p < nprocs; p++) {
            /* Read kloc[p] elements from file */
            char* buf = (p == 0) ? (char*) ap       // read directly into root process' output buffer
                                 : (char*) tmp;     // read into temporary buffer before sending to another process

            if(packed) {
                /* Packed storage, where the n*(n+1)/2 unique matrix elements are
                 * arranged in a flat array, in column-major upper-diagonal packed
                 * storage.  These elements will be read directly from file and
                 * distributed evenly among processes. */

                size_t nread = f.read(buf, kloc[p]*sizeof(real));
                if(nread != kloc[p]*sizeof(real)) {
                    opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd)\n", f.get_filename().c_str(), p, kloc[p], nread);
                    opsec_abort(1);
                }
            }
            else {
                /* Full matrix storage, where all n*n matrix elements are stored.
                 * Since the matrix is symmetric, we only need n*(n+1)/2 of them
                 * for subsequent calculations, so we discard the rest while
                 * reading from file. */

                int nsofar = 0;         // number of elements read so far, out of kloc[p]
                while(nsofar < kloc[p]) {
                    int nneeded = kloc[p] - nsofar;
                    int nleftincol = (col + 1 - row);
                    if(nneeded > nleftincol) {
                        /* Read in all elements left in column, increment column, and reset row */
                        size_t nread = f.read(buf, nleftincol*sizeof(real));
                        if(nread != nleftincol*sizeof(real)) {
                            opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd,%d,%d,%d)\n", f.get_filename().c_str(), p, kloc[p], nread, col, row, nleftincol);
                            opsec_abort(1);
                        }
                        buf += nleftincol*sizeof(real);
                        col += 1;
                        row = 0;
                    }
                    else {
                        /* Read only the elements that we need, and update row */
                        size_t nread = f.read(buf, nneeded*sizeof(real));
                        if(nread != nneeded*sizeof(real)) {
                            opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd,%d,%d,%d)\n", f.get_filename().c_str(), p, kloc[p], nread, col, row, nneeded);
                            opsec_abort(1);
                        }
                        row += nneeded;
                    }
                }
            }
        }

#if 0
        if(packed) {
            /* Packed storage, where the n*(n+1)/2 unique matrix elements are
             * arranged in a flat array, in column-major upper-diagonal packed
             * storage.  These elements will be read directly from file and
             * distributed evenly among processes. */

            for(int p = 0; p < nprocs; p++) {
                char* buf = (p == 0) ? (char*) ap       // read directly into root process' output buffer
                                     : (char*) tmp;     // read into temporary buffer before sending to another process
                size_t nread = f.read(buf, kloc[p]*sizeof(real));
                if(nread != kloc[p]*sizeof(real)) {
                    opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd)\n", f.get_filename().c_str(), p, kloc[p], nread);
                    opsec_abort(1);
                }
                if(p != 0) {
                    /* Send elements to process p */
                    MPI_Send(buf, kloc[p], REAL_MPI_TYPE, p, 0, MPI_COMM_WORLD);
                }
            }
        }
        else {
            /* Full matrix storage, where all n*n matrix elements are stored.
             * Since the matrix is symmetric, we only need n*(n+1)/2 of them
             * for subsequent calculations, so we discard the rest while
             * reading from file. */

            int col = 0, row = 0;       // current column and row number within file
            for(int p = 0; p < nprocs; p++) {
                char* buf = (p == 0) ? (char*) ap       // read directly into root process' output buffer
                                     : (char*) tmp;     // read into temporary buffer before sending to another process
                int nsofar = 0;         // number of elements read so far, out of kloc[p]
                while(nsofar < kloc[p]) {
                    int nneeded = kloc[p] - nsofar;
                    int nleftincol = (col + 1 - row);
                    if(nneeded > nleftincol) {
                        /* Read in all elements left in column, increment column, and reset row */
                        size_t nread = f.read(buf, nleftincol*sizeof(real));
                        if(nread != nleftincol*sizeof(real)) {
                            opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd,%d,%d,%d)\n", f.get_filename().c_str(), p, kloc[p], nread, col, row, nleftincol);
                            opsec_abort(1);
                        }
                        buf += nleftincol*sizeof(real);
                        col += 1;
                        row = 0;
                    }
                    else {
                        /* Read only the elements that we need, and update row */
                        size_t nread = f.read(buf, nneeded*sizeof(real));
                        if(nread != nneeded*sizeof(real)) {
                            opsec_error("ReadSymmetricMatrix: read error (%s,%d,%d,%zd,%d,%d,%d)\n", f.get_filename().c_str(), p, kloc[p], nread, col, row, nneeded);
                            opsec_abort(1);
                        }
                        row += nneeded;
                    }
                }
            }
        }
#endif
        free(tmp);
    }

    if(me != 0) {
        /* Receive elements from root process */
        MPI_Status status;
        MPI_Recv(ap, kloc[me], REAL_MPI_TYPE, 0, 0, MPI_COMM_WORLD, &status);
    }
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
    if(!cfg_has_keys(cfg, "estfile,evalfile,covfile,pixfile,Nmodes,packed_matrices", ",")) {
        fprintf(stderr, "estimate: missing configuration options\n");
        fputs(usage, stderr);
        opsec_exit(1);
    }

    const char* estfile = cfg_get(cfg, "estfile");
    const char* evalfile = cfg_get(cfg, "evalfile");
    const char* covfile = cfg_get(cfg, "covfile");
    const char* pixfile = cfg_get(cfg, "pixfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");
    int packed = cfg_get_int(cfg, "packed_matrices");

    std::string mixing = "inverse";
    if(cfg_has_key(cfg, "mixing"))
        mixing = cfg_get(cfg, "mixing");

    /* Load model, to get prior parameter values */
    Model* model = InitializeModel(cfg);
    int Nparams = model->NumParams();
    std::vector<real> param(Nparams);
    for(int n = 0; n < Nparams; n++)
        param[n] = model->GetParam(n);


    /* Number of coefficients needed for a symmetrix Nmodes-by-Nmodes matrix */
    int Ncoeff = Nmodes*(Nmodes+1)/2;

    /* Determine the local size and offset for each process */
    std::vector<int> locsize(nprocs), locdisp(nprocs);
    for(int p = 0; p < nprocs; p++) {
        locsize[p] = (Ncoeff/nprocs) + (p < (Ncoeff % nprocs));
        locdisp[p] = p*(Ncoeff/nprocs) + std::min(p, Ncoeff % nprocs);
    }
    int mylocsize = locsize[me];
    int mylocdisp = locdisp[me];


    /***** Memory allocation **************************************************/

    /* Signal-to-noise eigenvalues $\lambda_i$ will be distributed to all
     * processes. */
    real* lvalues = (real*) opsec_malloc(Nmodes*sizeof(real));

    /* Covariance matrix derivatives will be distributed across processes, with
     * process p getting locsize[p] of the Nmodes*(Nmodes+1)/2 elements.  The
     * kth value stored by process p (0 <= k < locsize[p]) will be element
     * number
     *   r = locdisp[p] + k
     * overall, corresponding to the matrix element D_{ij} with
     *   r = i*(i+1)/2 + j */
    real* dvalues = (real*) opsec_malloc(Nparams*mylocsize*sizeof(real));
    /* Fill with -1 to detect errors more easily */
    for(int m = 0; m < Nparams; ++m)
        for(int k = 0; k < mylocsize; ++k)
            dvalues[m*mylocsize + k] = -1;

    /* Residue matrix, symmetric, distributed over processes as above.
     * Initialized to zero. */
    real* rvalues = (real*) opsec_calloc(mylocsize, sizeof(real));

    /* Fisher matrix F_{mn}, quadratic estimates q_n, and bias component f_n. */
    real* fvalues = (real*) opsec_calloc(Nparams*(Nparams+2), sizeof(real));

    /* Starting and ending local indices for packed symmetric matrix */
    PackedMatrixIndex kbegin(Nmodes, mylocdisp);
    PackedMatrixIndex kend(Nmodes, mylocdisp + mylocsize);


    /***** Read data from file ************************************************/

    /* Read signal-to-noise eigenvalues $\lambda_i$ into lvalues array */
    SplitFile feval(evalfile, "r");
    ReadEvals(feval, Nmodes, lvalues);
    MyVector<real> lambda(Nmodes, lvalues);
    feval.close();

    /* Read covariance matrix derivatives $D_m$ into packed matrices */
    SplitFile fcov;
    if(me == 0) {
        /* Open file */
        fcov.open(covfile, "r");
        if(!fcov) {
            fprintf(stderr, "estimate-mpi: cannot open covfile (%s)\n", covfile);
            opsec_abort(1);
        }
        /* Read header */
        ReadCovHeader(fcov, Nparams, Nmodes, packed);
    }
    for(int m = 0; m < Nparams; m++)
        ReadSymmetricMatrix(fcov, Nmodes, &dvalues[m*mylocsize], packed);
    if(me == 0)
        fcov.close();

    /* Load matrix values into vectors */
    std::vector<MyVector<real> > D(Nparams);
    for(int m = 0; m < Nparams; m++)
        D[m] = MyVector<real>(mylocsize, &dvalues[m*mylocsize]);

    /* Print $C,m$ */
//    for(int m = 0; m < Nparams; m++)
//        for(PackedMatrixIndex k = kbegin; k < kend; ++k)
//            printf("Process %d: D[%d](%d,%d) = %g\n", me, m, k.i, k.j, D[m][k]);

    /* Read pixel values $y_i$ into yvalues */
    int Npixels;
    real* yvalues = ReadPixels(pixfile, &Npixels);
    if(yvalues == NULL || Npixels != Nmodes) {
        fprintf(stderr, "estimate: error reading pixel values\n");
        opsec_exit(1);
    }
    MyVector<real> y(Npixels, yvalues);

    /* Compute diagonal elements of covariance matrix,
     *   C_{ii} = 1 + \lambda_i */
    MyVector<real> C(Nmodes, lvalues);
    for(int i = 0; i < Nmodes; i++)
        C[i] = 1. + lambda[i];

    /* Compute residue matrix $R_{ij} = C_{ij} - \sum_n p_n C_{ij,n}$ */
    MyVector<real> R(mylocsize, rvalues);       // rvalues initialized to zero, so R_{ij} = 0
    for(PackedMatrixIndex k = kbegin; k < kend; ++k) {
        int i = k.i, j = k.j;
        if(i == j)
            R[k] += C[i];                   // R_{ij} += C_{ij}
        for(int n = 0; n < Nparams; n++)
            R[k] -= param[n] * D[n][k];     // R_{ij} -= p_n C_{ij,n}
    }

    /* Invert covariance matrix (it's diagonal and positive definite, so this is easy) */
    MyVector<real> Cinv(Nmodes, lvalues);
    for(int i = 0; i < Nmodes; i++)
        Cinv[i] = 1/C[i];

    /* Compute Fisher matrix $F_{mn} = (1/2) \Tr[C^{-1} C_{,m} C^{-1} C_{,n}]$ */
    MyMatrix<real> F(Nparams, Nparams, fvalues);
    #pragma omp parallel for
    for(int m = 0; m < Nparams; m++) {
        for(int n = 0; n <= m; n++) {
            real Fmn = 0;
            for(PackedMatrixIndex k = kbegin; k < kend; ++k) {
                int i = k.i, j = k.j;
                real mult = (i == j) ? 0.5 : 1.;      // multiplicative factor, 1 or 1/2
                Fmn += mult * Cinv[i] * D[m][k] * Cinv[j] * D[n][k];
            }
            F(m,n) = Fmn;
        }
    }
    for(int m = 0; m < Nparams; m++)
        for(int n = m+1; n < Nparams; n++)
            F(m,n) = F(n,m);

    /* Compute bias corrections $f_n = (1/2) Tr[C^{-1} C_{,n} C^{-1} R]$ */
    MyVector<real> f(Nparams, &fvalues[Nparams*Nparams]);
    #pragma omp parallel for
    for(int n = 0; n < Nparams; n++) {
        real fn = 0;
        for(PackedMatrixIndex k = kbegin; k < kend; ++k) {
            int i = k.i, j = k.j;
            real mult = (i == j) ? 0.5 : 1.;      // multiplicative factor, 1 or 1/2
            fn += mult * Cinv[i] * D[n][k] * Cinv[j] * R[k];
        }
        f[n] = fn;
    }

    /* Compute quadratic estimators $q_n = (1/2) y^T C^{-1} C_{,n} C^{-1} y$ */
    MyVector<real> q(Nparams, &fvalues[Nparams*(Nparams+1)]);
    #pragma omp parallel for
    for(int n = 0; n < Nparams; n++) {
        real qn = 0;
        for(PackedMatrixIndex k = kbegin; k < kend; ++k) {
            int i = k.i, j = k.j;
            real mult = (i == j) ? 0.5 : 1.;      // multiplicative factor, 1 or 1/2
            qn += mult * y[i] * Cinv[i] * D[n][k] * Cinv[j] * y[j];
        }
        q[n] = qn;
    }

#ifdef OPSEC_USE_MPI
    /* Reduce all results to root process */
    if(me == 0)
        MPI_Reduce(MPI_IN_PLACE, fvalues, Nparams*(Nparams+2), REAL_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
        MPI_Reduce(fvalues, NULL, Nparams*(Nparams+2), REAL_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

    /* Free storage that is no longer needed */
    free(lvalues);
    free(yvalues);
    free(dvalues);
    free(rvalues);
    if(me != 0)
        free(fvalues);

    /* Let the root process finish the rest */
    if(me == 0) {
        /* Memory allocation for matrices $M$, $W$, etc. */
        real* mvalues = (real*) opsec_malloc(Nparams*Nparams*sizeof(real));
        real* wvalues = (real*) opsec_malloc(Nparams*Nparams*sizeof(real));
        real* tvalues = (real*) opsec_malloc(Nparams*Nparams*sizeof(real));
        real* pvalues = (real*) opsec_malloc(Nparams*(Nparams+1)*sizeof(real));

        /* Choose mixing matrix */
        MyMatrix<real> M(Nparams, Nparams, mvalues);
        if(mixing.substr(0, 7) == "inverse") {
            inverse(F, M);              // M = F^{-1}
            if(mixing == "inverse sqrt")
                sqrtm(M, M);            // M = F^{-1/2}
        }
        else if(mixing == "identity") {
            for(int m = 0; m < Nparams; m++)
                for(int n = 0; n < Nparams; n++)
                    M(m,n) = (m == n);  // M = identity
        }
        else {
            fprintf(stderr, "estimate-mpi: unrecognized option: mixing = '%s'\n", mixing.c_str());
            opsec_exit(1);
        }

        /* Compute window matrix */
        MyMatrix<real> W(Nparams, Nparams, wvalues);
        multiply(M, F, W);      // W = M*F

        /* ... and normalize M so that the rows of W sum to unity */
        real rowsum;
        for(int m = 0; m < Nparams; m++) {
            rowsum = 0;
            for(int n = 0; n < Nparams; n++)
                rowsum += W(m,n);
            for(int n = 0; n < Nparams; n++)
                M(m,n) /= rowsum;
        }
        multiply(M, F, W);      // W = M*F

        /* Compute power estimates */
        MyVector<real> phat(Nparams, pvalues);
        axpy(-1., f, q);         // q = q - f
        multiply(M, q, phat);   // \hat{p} = M*(q-f)

        /* ... and their covariance matrix */
        MyMatrix<real> tmp(Nparams, Nparams, tvalues);
        MyMatrix<real> pcov(Nparams, Nparams, &pvalues[Nparams]);
        multiply(M, F, tmp);
        multiply(tmp, M.transpose(), pcov);

        /* Print them out */
        printf("p_n          error\n");
        for(int n = 0; n < Nparams; n++)
            printf("%e %e\n", phat[n], sqrt(pcov(n,n)));

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
        time_t now;
        time(&now);
        fprintf(fest, "# Generated %s\n", ctime(&now));

        Config opts = cfg_new();
        cfg_set_int(opts, "Nparams", Nparams);

        abn_write(fest, &phat[0], Nparams, REAL_FMT, opts);
        abn_write(fest, &pcov(0,0), Nparams*Nparams, REAL_FMT, NULL);
        abn_write(fest, &F(0,0), Nparams*Nparams, REAL_FMT, NULL);
        abn_write(fest, &W(0,0), Nparams*Nparams, REAL_FMT, NULL);
        abn_write(fest, &M(0,0), Nparams*Nparams, REAL_FMT, NULL);
        abn_write(fest, &q[0], Nparams, REAL_FMT, NULL);
        abn_write(fest, &f[0], Nparams, REAL_FMT, NULL);

//        F.print();
//        W.print();
//        M.print();

        cfg_destroy(opts);
        fclose(fest);
        free(mvalues);
        free(wvalues);
        free(pvalues);
        free(tvalues);
        free(fvalues);
    }

    /* Clean up and exit */
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif
    return 0;
}
