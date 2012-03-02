/* klt
 *
 * Compute the signal covariance matrix in the cell function basis, then
 * compute its eigenvector and eigenvalues.  The eigenvectors define the KL
 * mode functions, while the eigenvalues define the (diagonal elements of the)
 * signal matrix in the KL basis. */

/* The signal matrix is computed first in the cell basis, distributed among
 * processes as illustrated below:
 *
 *                           n  
 *             +----------------------------+
 * Process 0:  |                            |  nloc
 *             +----------------------------+
 * Process 1:  |                            |
 *             +----------------------------+
 *             :                            :
 *             :                            :
 *             +----------------------------+
 * nprocs-1    |                            |
 *             +----------------------------+
 *
 * Each nloc-by-n block of the matrix is stored locally as an array S in
 * row-major (i.e. C-style) order.  That is, each of the blocks above would be
 * filled in with array values as
 *     S[0]           S[1]         S[2]    ...  S[n-1]
 *     S[n]           S[n+1]       S[n+1]  ...  S[2*n-1]
 *      :               :            :     ...     :
 *     S[(nloc-1)*n]  ...                  ...  S[nloc*n-1]
 *
 * Once the local blocks of the signal matrix are computed, PARPACK is invoked
 * to compute the largest Nmodes eigenvalues and eigenvectors.
 */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <string>
#include <vector>
using std::string;
using std::vector;

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "Cell.h"
#include "Model.h"
#include "Survey.h"
#include "XiFunc.h"
#include "abn.h"
#include "cfg.h"
#include "eig.h"
#include "opsec.h"
#include "sig.h"

const char* usage =
    "Usage: klt [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h             Display this help\n"
    "  -c FILE        Read additional configuration options from FILE\n"
    "Configuration options:\n"
    "  cellfile=FILE  Read list of cells from FILE (required)\n"
    "  Nmodes=NUM     Compute NUM mode functions of greatest S/N (required)\n"
    "  modefile=FILE  Write mode function coefficients to FILE (required)\n"
    "  evalfile=FILE  Write S/N eigenvalues to FILE (required)\n"
    "  model=NAME     Use specified model (required)\n"
    "  coordsys=TYPE  Coordinate system; either spherical or cartesian (required)\n"
    "  Nr,Nmu,Nphi    Spherical grid (required if coordsys = spherical)\n"
    "  Nx,Ny,Nz       Cartesian grid (required if coordsys = cartesian)\n";

int main(int argc, char* argv[]) {
    /* Initialize MPI, if available */
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    Config cfg = opsec_init(argc, argv, usage);

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "coordsys,cellfile,Nmodes,modefile,evalfile", ",")) {
        opsec_error("klt: missing configuration options\n");
        if(me == 0) { fputs(usage, stderr); fflush(stderr); }
        opsec_exit(1);
    }

    int coordsys = cfg_get_enum(cfg, "coordsys", "cartesian", CoordSysCartesian,
                                                 "spherical", CoordSysSpherical,
                                                 "", -1);
    if(coordsys == -1) {
        opsec_error("klt: missing or invalid config option: coordsys = %s\n", cfg_get(cfg, "coordsys"));
        opsec_exit(1);
    }

    const char* cellfile = cfg_get(cfg, "cellfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");
    const char* modefile = cfg_get(cfg, "modefile");
    const char* evalfile = cfg_get(cfg, "evalfile");

    /* Load model */
    Model* model = InitializeModel(cfg);
    if(!model)
        opsec_exit(1);
    XiFunc xi = model->GetXi();

    /* Load survey */
    Survey* survey = InitializeSurvey(cfg);
    if(!survey)
        opsec_exit(1);

    /* Open output files on root process */
    FILE* fmode = NULL;
    FILE* feval = NULL;
    if(me == 0) {
        fmode = fopen(modefile, "w");
        if(fmode == NULL) {
            opsec_error("klt: could not open '%s' for writing\n", modefile);
            opsec_exit(1);
        }
        feval = fopen(evalfile, "w");
        if(feval == NULL) {
            opsec_error("klt: could not open '%s' for writing\n", evalfile);
            opsec_exit(1);
        }
    }

    /* Read cells from file */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);

    /* Determine the local problem lengths for each process */
    vector<int> locsizes(nprocs), locdisps(nprocs);
    for(int i = 0; i < nprocs; i++) {
        locsizes[i] = (Ncells/nprocs) + (i < (Ncells % nprocs));
        locdisps[i] = (i == 0) ? 0 : locdisps[i-1] + locsizes[i-1];
    }
    assert(locdisps[nprocs-1] + locsizes[nprocs-1] == Ncells);

    /* Total number of cells */
    int n = Ncells;

    /* Number of rows of n-by-n matrix that this process is responsible for */
    int nloc = locsizes[me];

    /* Index of first row that this process is responsible for */
    int amin = locdisps[me];

    /* Memory for local nloc-by-n block of signal matrix */
    real* S = NULL;

    opsec_info("Computing signal matrix elements in cell basis...\n");

    MPI_Barrier(MPI_COMM_WORLD);

    if(coordsys == CoordSysSpherical) {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi", ",")) {
            opsec_error("klt: must provide config options N{r,mu,phi}\n");
            opsec_exit(1);
        }
        int Nr = cfg_get_int(cfg, "Nr");
        int Nmu = cfg_get_int(cfg, "Nmu");
        int Nphi = cfg_get_int(cfg, "Nphi");
        S = ComputeSignalMatrixS(n, nloc, amin, cells, Nr, Nmu, Nphi, xi, survey);
    }
    else if(coordsys == CoordSysCartesian) {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz", ",")) {
            opsec_error("klt: must provide config options N{x,y,z}\n");
            opsec_exit(1);
        }
        int Nx = cfg_get_int(cfg, "Nx");
        int Ny = cfg_get_int(cfg, "Ny");
        int Nz = cfg_get_int(cfg, "Nz");
        S = ComputeSignalMatrixC(n, nloc, amin, cells, Nx, Ny, Nz, xi, survey);
    }

    /* Number of eigenvalues to compute */
    int nev = std::min(Nmodes, n-1);

    /* Number of working basis vectors to use by ARPACK */
    int ncv = std::min(2*nev, n);

    /* Local memory for eigenvalues and eigenvectors */
    real* evals = (real*) opsec_malloc(nev*sizeof(real));
    real* modes = (real*) opsec_malloc(nloc*ncv*sizeof(real));

    opsec_info("Computing eigenmodes of signal matrix...\n");

#ifdef OPSEC_USE_MPI
    int nconv = peig(MPI_COMM_WORLD, n, nloc, nev, ncv, S, evals, modes);
#else
    int nconv = eig(n, nev, ncv, S, evals, modes);
#endif

    /* Free signal matrix memory as soon as possible */
    free(S);

    opsec_info("Writing eigenmodes and eigenvalues to file...\n");

    /* Prepare headers in output files */
    real* v = NULL;     // memory for a full eigenvector on root process
    if(me == 0) {

        v = (real*) opsec_malloc(n*sizeof(real));

        Config opts = cfg_new();
        cfg_set(opts, "cellfile", cellfile);
        cfg_set_int(opts, "Ncells", Ncells);
        cfg_set_int(opts, "Nmodes", Nmodes);

        fprintf(feval, "# Signal matrix eigenvalues\n");
        abn_write_header(feval, nev, REAL_FMT, opts);

        fprintf(fmode, "# Signal matrix eigenmodes\n");
        fprintf(fmode, "# Each mode is stored consecutively, as below:\n");
        fprintf(fmode, "#   v_1,1   v_1,2   ... v_1,n\n");
        fprintf(fmode, "#   ...\n");
        fprintf(fmode, "#   v_nev,1 v_nev,2 ... v_nev,n\n");
        abn_write_header(fmode, nev*n, REAL_FMT, opts);

        cfg_destroy(opts);
    }

    for(int i = 1; i <= nconv; i++) {
        /* (The strange indexing below comes from the fact that, while PARPACK
         * returns the 'nev' largest eigenvalues, it stores this list in
         * _reverse_ order, i.e. smallest to largest.  We want to write the
         * spectrum to file with the largest eigenvalue first.) */
#ifdef OPSEC_USE_MPI
        MPI_Gatherv(modes + (nev-i)*nloc, nloc, REAL_MPI_TYPE, v, &locsizes[0], &locdisps[0], REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
#else
        memcpy(v, modes + (nev-i)*n, n*sizeof(real));
#endif
        if(me == 0) {
            fwrite(&evals[nev-i], sizeof(real), 1, feval);
            fwrite(v, sizeof(real), n, fmode);
        }
    }

    /* Clean up nicely */
    if(me == 0) {
        fclose(feval);
        fclose(fmode);
        free(v);
    }
    free(evals);
    free(modes);
    free(cells);
    delete model;
    delete survey;
    cfg_destroy(cfg);
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
