/* comma
 *
 * Compute the derivatives of the covariance matrix with respect to the model
 * parameters.  First compute this matrix in the cell function basis, then
 * use the KL mode coefficients to project onto the KL basis.
 *
 * Since this code currently needs to store an entire Nmodes-by-Nmodes matrix
 * in memory for each process, it's probably best to run it with only one or
 * two processes per node.  If the BLAS implementation supports OpenMP-level
 * parallelism (which MKL and ACML both claim to do), the extra cores will be
 * made use of in the calls to gemv and gemm. */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <string>
#include <vector>
using std::string;
using std::vector;

#ifdef OPSEC_USE_MPI
# include <mpi.h>
#endif

#include "Cell.h"
#include "Model.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "sig.h"

const char* usage =
    "Usage: comma [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h             Display this help\n"
    "  -c FILE        Read additional configuration options from FILE\n"
    "Configuration options:\n"
    "  cellfile=FILE  Read list of cells from FILE (required)\n"
    "  modefile=FILE  Read KL modes from FILE (required)\n"
    "  covfile=FILE   Write covariance matrices to FILE (required)\n"
    "  Nmodes=NUM     Number of KL modes\n"
    "  model=NAME     Use specified model (required)\n"
    "  coordsys=TYPE  Coordinate system; either spherical or cartesian (required)\n"
    "  Nr,Nmu,Nphi    Spherical grid (required if coordsys = spherical)\n"
    "  Nx,Ny,Nz       Cartesian grid (required if coordsys = cartesian)\n";

/* Project (local) signal matrix from cell basis to KL basis.  S is the signal
 * matrix, nloc-by-Ncells stored in C (row-major) order.  B is the projection
 * matrix, Nmodes-by-Ncells also stored in C order.  out is the resulting
 * matrix, Nmodes-by-Nmodes stored in C order. */
void LocalProject(int Ncells, int nloc, int Nmodes, int amin, real* B, real* S, real* out) {
    char trans = 'T', notrans = 'N';
    real alpha = 1.0, beta = 0.0;

    /* Compute product tmp = S * B^t */
    real* tmp = (real*) malloc(nloc*Nmodes*sizeof(real));
    blas_gemm(&trans, &notrans, &nloc, &Nmodes, &Ncells, &alpha, S, &Ncells, B, &Ncells, &beta, tmp, &nloc);

    /* Compute product out = Bloc * tmp, but in C order rather than Fortran order.
     * - B is Nmodes-by-Ncells stored in C order, so Ncells-by-Nmodes in Fortran order
     * - tmp is nloc-by-Nmodes stored in Fortran order, so Nmodes-by-nloc in C order
     * - Computing out = trans(tmp)*B will give the right answer in C order (trust me) */
    real* Bloc = B + amin;
    blas_gemm(&trans, &notrans, &Nmodes, &Nmodes, &nloc, &alpha, tmp, &nloc, Bloc, &Ncells, &beta, out, &Nmodes);

    free(tmp);
}

int main(int argc, char* argv[]) {
    int Ncells, Nmodes, Nparams;
    const char *cellfile, *modefile, *covfile;
    FILE* fout = NULL;

    /* Initialize MPI */
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);
#endif

    Config cfg = opsec_init(argc, argv, usage);

    /* Debugging... */
//    if(me == 0) {
//        printf("# Config options\n");
//        cfg_write(cfg, stdout);
//    }

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "coordsys,cellfile,modefile,covfile,Nmodes", ",")) {
        if(me == 0) fprintf(stderr, "comma: missing configuration options\n");
        if(me == 0) fputs(usage, stderr);
        opsec_exit(1);
    }
    cellfile = cfg_get(cfg, "cellfile");
    modefile = cfg_get(cfg, "modefile");
    covfile = cfg_get(cfg, "covfile");
    Nmodes = cfg_get_int(cfg, "Nmodes");

    /* Determine output format; default to single file */
    string comma_output = cfg_has_key(cfg, "comma_output") ? cfg_get(cfg, "comma_output") : "single";
    bool multifile = (comma_output.substr(0,5) == "multi");

#if 0
    /* Open output file first to make sure we'll be able to write */
    if(me == 0) {
        fout = fopen(covfile, "w");
        if(fout == NULL) {
            fprintf(stderr, "comma: could not open '%s' for writing\n", covfile);
            exit_comma(1);
        }
    }
#endif

    /* Decide on coordinate system */
    string coordsys(cfg_get(cfg, "coordsys"));
    if(coordsys != "spherical" && coordsys != "cartesian") {
        fprintf(stderr, "comma: missing or invalid config option: coordsys = %s\n", coordsys.c_str());
        opsec_exit(1);
    }

    int Nr, Nmu, Nphi, Nx, Ny, Nz;
    if(coordsys == "spherical") {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi", ",")) {
            if(me == 0)
                fprintf(stderr, "comma: need Nr,Nmu,Nphi for spherical coordinates\n");
            opsec_exit(1);
        }
        Nr = cfg_get_int(cfg, "Nr");
        Nmu = cfg_get_int(cfg, "Nmu");
        Nphi = cfg_get_int(cfg, "Nphi");
    }
    else if(coordsys == "cartesian") {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz", ",")) {
            if(me == 0)
                fprintf(stderr, "comma: need Nx,Ny,Nz for cartesian coordinates\n");
            opsec_exit(1);
        }
        Nx = cfg_get_int(cfg, "Nx");
        Ny = cfg_get_int(cfg, "Ny");
        Nz = cfg_get_int(cfg, "Nz");
    }

    /* Get model */
    Model* model = InitializeModel(cfg);
    if(!model)
        opsec_exit(1);
    Nparams = model->NumParams();

    /* Read cells from file */
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);

#if 0
    /* Prepare modes file on process 0 */
    FILE* fmodes = NULL;
    long first_mode = 0;
    if(me == 0) {
        int Nmodes_check, Ncells_check;
        fmodes = ReadModesHeader(modefile, &Nmodes_check, &Ncells_check);
        first_mode = ftell(fmodes);
        if(Nmodes != Nmodes_check || Ncells != Ncells_check) {
            fprintf(stderr, "comma: inconsistent modes file\n");
            exit_comma(1);
        }
    }
#endif

    int Nmodes_check, Ncells_check;
    real* B = ReadModes(modefile, &Nmodes_check, &Ncells_check);
    if(Nmodes != Nmodes_check || Ncells != Ncells_check) {
        fprintf(stderr, "comma: inconsistent modes file\n");
        opsec_exit(1);
    }

    /* Determine the local problem lengths for each process */
    vector<int> locsizes(nprocs), locdisps(nprocs);
    for(int i = 0; i < nprocs; i++) {
        locsizes[i] = (Ncells/nprocs) + (i < (Ncells % nprocs));
        locdisps[i] = (i == 0) ? 0 : locdisps[i-1] + locsizes[i-1];
    }
    assert(locdisps[nprocs-1] + locsizes[nprocs-1] == Ncells);

    int nloc = locsizes[me];
    int amin = locdisps[me];

    /* Each process computes partial contributions to the matrix elements
     *   \tilde{C}_{ij,m} = \sum_{a,b} B_{ia} B_{jb} C_{ab,m}
     * for 0 <= i,j < Nmodes and 0 <= m < Nparams.  Process p will compute
     * C_{ab,m} for 0 <= m < Nparams, amin <= a < amin + nloc, and 0 <= b <
     * Ncells, and perform the appropriate partial sum.  The partial results
     * are then gathered and finally summed by the root process, which writes
     * the full matrices out to file. */


    /* Initialize correlation functions */
    vector<XiFunc> xip(Nparams);        // xip[m] = $\xi_{,m}$
    for(int m = 0; m < Nparams; m++)
        xip[m] = model->GetXiDeriv(m);

    /* For single file output, prepare header */
    if(!multifile && me == 0) {
        fout = fopen(covfile, "wb");
        if(fout == NULL) {
            fprintf(stderr, "comma: could not open '%s' for writing\n", covfile);
            opsec_exit(1);
        }

        fprintf(fout, "# Derivatives of the covariance matrix in KL basis, C_{ij,m}.\n");
        fprintf(fout, "# Single file output\n");
        fprintf(fout, "# Binary data consists of %d consecutive %dx%d matrices.\n", Nparams, Nmodes, Nmodes);

        Config opts = cfg_new();
        cfg_set_int(opts, "Nmodes", Nmodes);
        cfg_set_int(opts, "Nparams", Nparams);
        abn_write_header(fout, Nparams*Nmodes*Nmodes, "d", opts);
        cfg_destroy(opts);
    }

    if(me == 0)
        printf("Computing parameter derivatives of covariance matrix...\n");


#if 0
    /* Largest instantaneous memory allocation is either
     *   (2*nloc*Nmodes + nloc*Ncells)*sizeof(real) bytes
     * or
     *   (2*nloc*Nmodes + Nmodes*Nmodes)*sizeof(real) bytes,
     * whichever is larger. */
#endif

    /* Compute matrix derivatives one parameter at a time */
    for(int m = 0; m < Nparams; m++) {
        if(me == 0) {
            printf(" Parameter %d/%d\n", m+1, Nparams);
            fflush(stdout);
        }

        /* For multi-file output, prepare mth output file */
        if(multifile && me == 0) {
            char filename[512];
            snprintf(filename, 512, "%s.%03d", covfile, m);
            fout = fopen(filename, "wb");
            if(fout == NULL) {
                fprintf(stderr, "comma: could not open '%s' for writing\n", filename);
                opsec_exit(1);
            }

            fprintf(fout, "# Derivatives of the covariance matrix in KL basis, C_{ij,m}.\n");
            fprintf(fout, "# Multi-file output; parameter number %d.\n", m);
            fprintf(fout, "# Binary data consists of a single %dx%d matrix.\n", Nmodes, Nmodes);

            Config opts = cfg_new();
            cfg_set_int(opts, "m", m);
            cfg_set_int(opts, "Nmodes", Nmodes);
            cfg_set_int(opts, "Nparams", Nparams);
            abn_write_header(fout, Nmodes*Nmodes, "d", opts);
            cfg_destroy(opts);
        }

#if 0
        /* Set file pointer to the first KL eigenmode */
        if(me == 0)
            fseek(fmodes, first_mode, SEEK_SET);
#endif

        /* Compute local block of the matrix $C_{ab,m}$ in cell basis */
        real* S = NULL;
        if(coordsys == "spherical")
            S = ComputeSignalMatrixS(Ncells, nloc, amin, cells, Nr, Nmu, Nphi, xip[m]);
        else if(coordsys == "cartesian")
            S = ComputeSignalMatrixC(Ncells, nloc, amin, cells, Nx, Ny, Nz, xip[m]);

        if(me == 0) {
            printf("  Projecting onto KL basis...\n");
            fflush(stdout);
        }
        real* Sp = (real*)opsec_malloc(Nmodes*Nmodes*sizeof(real));
        LocalProject(Ncells, nloc, Nmodes, amin, B, S, Sp);

#if 0
        /* Storage space for a single KL eigenmode */
        real* mode = (real*)opsec_malloc(Ncells*sizeof(real));
        /* Storage for the local submatrix of the projection matrix B */
        real* B = (real*)opsec_malloc(Nmodes*nloc*sizeof(real));
        /* Storage for the temporary nloc-by-Nmodes submatrix of S.B^t */
        real* SBt = (real*)opsec_malloc(nloc*Nmodes*sizeof(real));

        /* Read KL eigenmodes from disk one mode at a time */
        for(int i = 0; i < Nmodes; i++) {
            if(me == 0) {
                size_t nread = fread(mode, sizeof(real), Ncells, fmodes);
                assert(nread == Ncells);
                printf("  Reading mode %d/%d\n", i+1, Nmodes);
            }
#ifdef OPSEC_USE_MPI
            MPI_Bcast(mode, Ncells, REAL_MPI_TYPE, 0, comm);
#endif

            /* Save the columns of B that we will need to comptue S' = B.S.B^t */
            for(int a = amin; a < amin + nloc; a++)
                B[nloc*i + (a-amin)] = mode[a];

            /* Compute local column of S.B^t using BLAS voodoo */
            static char trans = 'T';
            static real alpha = 1.0, beta = 0.0;
            static int incx = 1, incy = 1;
            blas_gemv(&trans, &Ncells, &nloc, &alpha, S, &Ncells, mode, &incx, &beta, &SBt[i], &Nmodes);
        }

        /* Free up some space */
        free(mode);
        free(S);

        /* Allocate Nmodes*Nmodes*sizeof(real) bytes of memory */
        real* Sp = (real*)opsec_malloc(Nmodes*Nmodes*sizeof(real));
        if(Sp == NULL) {
            fprintf(stderr, "comma: could not allocate %zd bytes for Sp on process %d\n", Nmodes*Nmodes*sizeof(real), me);
            exit_comma(1);
        }

        /* Compute S' = B.(S.B^t) using BLAS voodoo (hopefully parallelized by
         * BLAS implementation) */
        static char notrans = 'N';
        static real alpha = 1.0, beta = 0.0;
        blas_gemm(&notrans, &notrans, &Nmodes, &Nmodes, &nloc, &alpha, SBt, &Nmodes, B, &nloc, &beta, Sp, &Nmodes);

        free(B);
        free(SBt);
#endif

        if(me == 0) {
            printf("  Summing partial matrix contributions...\n");
            fflush(stdout);
        }

        /* Gather and sum contributions to $\tilde{C}_{ij,m}$ */
#ifdef OPSEC_USE_MPI
        if(me == 0)
            MPI_Reduce(MPI_IN_PLACE, Sp, Nmodes*Nmodes, REAL_MPI_TYPE, MPI_SUM, 0, comm);
        else
            MPI_Reduce(Sp, NULL, Nmodes*Nmodes, REAL_MPI_TYPE, MPI_SUM, 0, comm);
#endif

        /* Write matrix to file */
        if(me == 0)
            fwrite(Sp, sizeof(double), Nmodes*Nmodes, fout);

        free(S);
        free(Sp);

        /* For multi-file output, prepare for next file */
        if(multifile && me == 0) {
            fclose(fout);
            fout = NULL;
        }
    }

    /* Clean up nicely */
    if(!multifile && me == 0)
        fclose(fout);
    free(cells);
    delete model;
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
