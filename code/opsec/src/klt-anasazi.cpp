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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBasicSort.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#  include <Teuchos_RawMPITraits.hpp>
#endif

#include "Cell.h"
#include "Model.h"
#include "MyMultiVec.h"
#include "MySignalMatrixOperator.h"
#include "SplitFile.h"
#include "Survey.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "slp.h"

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

    /* Parse command line and configuration options */
    Config cfg = opsec_init(argc, argv, usage);

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "coordsys,cellfile,Nmodes,modefile,evalfile", ",")) {
        if(me == 0) opsec_error("klt: missing configuration options\n");
        if(me == 0) fputs(usage, stderr);
        opsec_exit(1);
    }

    int coordsys = cfg_get_enum(cfg, "coordsys", "cartesian", CoordSysCartesian,
                                                 "spherical", CoordSysSpherical,
                                                 "", -1);
    if(coordsys == -1) {
        opsec_error("klt: missing or invalid config option: coordsys = %s\n", cfg_get(cfg, "coordsys"));
        opsec_exit(1);
    }

    int N1, N2, N3;
    if(coordsys == CoordSysCartesian) {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz", ",")) {
            fprintf(stderr, "klt: must provide config options N{x,y,z}\n");
            opsec_exit(1);
        }
        N1 = cfg_get_int(cfg, "Nx");
        N2 = cfg_get_int(cfg, "Ny");
        N3 = cfg_get_int(cfg, "Nz");
    }
    else if(coordsys == CoordSysSpherical) {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi", ",")) {
            fprintf(stderr, "klt: must provide config options N{r,mu,phi}\n");
            opsec_exit(1);
        }
        N1 = cfg_get_int(cfg, "Nr");
        N2 = cfg_get_int(cfg, "Nmu");
        N3 = cfg_get_int(cfg, "Nphi");
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
    SplitFile fmode;
    SplitFile feval;
    if(me == 0) {
        fmode.open(modefile, "w");
        if(!fmode.isopen()) {
            opsec_error("klt: could not open '%s' for writing\n", modefile);
            opsec_abort(1);
        }
        feval.open(evalfile, "w");
        if(!feval.isopen()) {
            opsec_error("klt: could not open '%s' for writing\n", evalfile);
            opsec_abort(1);
        }
    }

    /* Read cells from file */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);

    if(Nmodes > Ncells) {
        if(me == 0) opsec_error("klt: requesting %d modes from a %d dimensional space\n", Nmodes, Ncells);
        opsec_exit(1);
    }

    /* Number of eigenvalues to compute */
    int nev = Nmodes;

    /* Block size to use for block Krylov-Schur method */
    int block_size = 1;

    /* Number of blocks to use */
    int num_blocks = std::min(2*nev/block_size, Ncells);

    /* Convergence tolerance */
    double tol = 1e-5;

    /* Maximum number of restarts */
    int max_restarts = 1000;

    /* Which eigenvalues to solve for (LM = largest magnitude) */
    std::string which = "LM";

    /* Verbosity level */
    int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary
                  + Anasazi::TimingDetails + Anasazi::StatusTestDetails;

    /* Construct parameter list for Eigensolver */
    Teuchos::ParameterList params;
    params.set("Verbosity", verbosity);
    params.set("Block Size", block_size);
    params.set("Num Blocks", num_blocks);
    params.set("Maximum Restarts", max_restarts);
    params.set("Convergence Tolerance", tol);
    params.set("Which", which);

    typedef real ST;
    typedef Anasazi::MultiVec<ST> MV;
    typedef Anasazi::Operator<ST> OP;
    typedef Anasazi::MultiVecTraits<ST,MV> MVT;
    typedef Anasazi::OperatorTraits<ST,MV,OP> OPT;

    if(me == 0)
        opsec_info("Initializing eigenproblem...\n");

    /* Construct initial vector */
    Teuchos::RCP<MyMultiVec<real> > ivec = Teuchos::rcp(new MyMultiVec<real>(Ncells, block_size));
    ivec->MvRandom();

    /* Construct operator */
    Teuchos::RCP<MySignalMatrixOperator<real> > op
        = Teuchos::rcp(new MySignalMatrixOperator<real>(coordsys, N1, N2, N3, Ncells, cells, xi, survey));

    /* Create the eigenproblem */
    Teuchos::RCP<Anasazi::BasicEigenproblem<ST, MV, OP> > problem =
        Teuchos::rcp(new Anasazi::BasicEigenproblem<ST, MV, OP>(op, ivec));
    problem->setHermitian(true);
    problem->setNEV(nev);
    bool problem_okay = problem->setProblem();
    if(!problem_okay && me == 0) {
        opsec_error("klt: error setting eigenproblem\n");
        opsec_abort(1);
    }

    if(me == 0)
        opsec_info("Initializing eigensolver...\n");

    /* Initialize the block-Arnoldi solver */
    Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> solmgr(problem, params);

    if(me == 0)
        opsec_info("Solving eigenproblem...\n");

    /* Solve the problem to the specified tolerances or length */
    Anasazi::ReturnType retcode = solmgr.solve();
    if(retcode != Anasazi::Converged && me == 0) {
        opsec_error("klt: SolverManager failed to converge\n");
    }

    if(me == 0)
        opsec_info("Reading solution...\n");

    /* Get the eigenvalues and eigenvectors from the eigenproblem */
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    std::vector<Anasazi::Value<ST> > evals = sol.Evals;
    Teuchos::RCP<MV> mvevecs = sol.Evecs;
    MyMultiVec<ST>* evecs = dynamic_cast<MyMultiVec<ST>*>(mvevecs.get());
    std::vector<int> index = sol.index;
    int nconv = sol.numVecs;

    if(nconv < nev) {
        if(me == 0)
            opsec_error("klt: not all eigenvalues converged: nev = %d, nconv = %d\n", nev, nconv);
        nev = nconv;
    }

    /* Prepare headers in output files */
    if(me == 0) {
        opsec_info("Writing eigenmodes and eigenvalues to file...\n");

        Config opts = cfg_new();
        cfg_set(opts, "cellfile", cellfile);
        cfg_set_int(opts, "Ncells", Ncells);
        cfg_set_int(opts, "Nmodes", Nmodes);

        /* Write evalfile header */
        fprintf(feval, "# Signal matrix eigenvalues\n");
        abn_write_header(feval, nev, REAL_FMT, opts);

        /* Write modefile header */
        fprintf(fmode, "# Signal matrix eigenmodes\n");
        fprintf(fmode, "# Each mode is stored consecutively, as below:\n");
        fprintf(fmode, "#   v_1,1   v_1,2   ... v_1,n\n");
        fprintf(fmode, "#   ...\n");
        fprintf(fmode, "#   v_nev,1 v_nev,2 ... v_nev,n\n");
        abn_write_header(fmode, nev*Ncells, REAL_FMT, opts);

        cfg_destroy(opts);
    }

    std::vector<real> v;        // memory for a full eigenvector on root process
    if(me == 0)
        v.resize(Ncells);

    /* Determine the local problem length and offset for each process */
    std::vector<int> locsizes(nprocs), locdisps(nprocs);
    for(int p = 0; p < nprocs; p++) {
        locsizes[p] = (Ncells/nprocs) + (p < (Ncells % nprocs));
        locdisps[p] = p*(Ncells/nprocs) + std::min(p, Ncells % nprocs);
    }

    for(int j = 0; j < nev; j++) {
        real* vloc = evecs->vector(j);

#ifdef OPSEC_USE_MPI
        MPI_Gatherv(vloc, locsizes[me], REAL_MPI_TYPE,
                    &v[0], &locsizes[0], &locdisps[0], REAL_MPI_TYPE,
                    0, MPI_COMM_WORLD);
#else
        memcpy(&v[0], vloc, Ncells*sizeof(real));
#endif
        if(me == 0) {
            feval.write((char*) &evals[j], sizeof(real));
            fmode.write((char*) &v[0], Ncells*sizeof(real));
        }
    }

    /* Clean up nicely */
    if(me == 0) {
        feval.close();
        fmode.close();
    }
    free(cells);
    delete model;
    delete survey;
    cfg_destroy(cfg);
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
