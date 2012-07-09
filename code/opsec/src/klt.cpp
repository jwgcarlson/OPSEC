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
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
using std::string;
using std::vector;

#include <mpi.h>

#include "AnasaziSolver.h"
#include "Cell.h"
#include "Model.h"
#include "ParpackSolver.h"
#include "SplitFile.h"
#include "Survey.h"
#include "XiFunc.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "sig.h"
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
    "  solver=NAME    Use specified eigensolver (currently parpack or anasazi)\n"
;

/* Helper class for computing local signal matrix elements. */
class SignalMatrixFactory : public MatrixFactory {
public:
    SignalMatrixFactory(const XiFunc& xi, Survey* survey, int Ncells, const Cell* cells, double epsrel = 1e-5, double epsabs = 1e-10)
        : xi(xi), survey(survey), Ncells(Ncells), cells(cells), epsrel(epsrel), epsabs(epsabs)
    {}

    void ComputeMatrixValues(const vector<int>& rows, const vector<int>& cols, real* values, int lld) {
        ComputeSignalMatrixBlock(Ncells, rows, cols, values, lld, xi, survey, cells, epsrel, epsabs);
    }

private:
    const XiFunc& xi;
    Survey* survey;
    int Ncells;
    const Cell* cells;
    double epsrel, epsabs;
};


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
    if(cfg_missing_keys(cfg, "cellfile,Nmodes,modefile,evalfile")) {
        if(me == 0) opsec_error("klt: missing configuration options\n");
        if(me == 0) fputs(usage, stderr);
        opsec_exit(1);
    }

    const char* cellfile = cfg_get(cfg, "cellfile");
    const char* modefile = cfg_get(cfg, "modefile");
    const char* evalfile = cfg_get(cfg, "evalfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");

    /* Read cells from file */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);

    if(Nmodes > Ncells) {
        if(me == 0) opsec_error("klt: requesting %d modes from a %d-dimensional space\n", Nmodes, Ncells);
        opsec_exit(1);
    }

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

    /* Create matrix factory */
    double epsrel = cfg_has_key(cfg, "sig.epsrel") ? cfg_get_double(cfg, "sig.epsrel") : 1e-5;
    double epsabs = cfg_has_key(cfg, "sig.epsabs") ? cfg_get_double(cfg, "sig.epsabs") : 1e-10;
    SignalMatrixFactory matfact(xi, survey, Ncells, cells, epsrel, epsabs);

    /* Select eigensolver */
    Solver* solver = NULL;
    string s = cfg_get(cfg, "solver");
    if(s == "parpack") {
#ifdef HAVE_PARPACK
        /* TODO: parse klt.parpack.* options and pass to solver */
        solver = new ParpackSolver(Ncells, matfact);
#else
        if(me == 0)
            opsec_error("klt: no support for PARPACK solver, try recompiling OPSEC\n");
        opsec_exit(1);
#endif
    }
    else if(s == "anasazi") {
#ifdef HAVE_ANASAZI
        /* TODO: parse klt.anasazi.* options and pass to solver */
        solver = new AnasaziSolver(Ncells, matfact);
#else
        if(me == 0)
            opsec_error("klt: no support for Anasazi solver, try recompiling OPSEC\n");
        opsec_exit(1);
#endif
    }
    else {
        if(me == 0)
            opsec_error("klt: unrecognized solver: %s\n", s.c_str());
        opsec_exit(1);
    }

    /* Solve for Nmodes eigenvalues/eigenvectors */
    int nconv = solver->Solve(Nmodes);

    /* Check for errors */
    if(nconv < Nmodes) {
        if(me == 0)
            opsec_error("klt: not all eigenvalues converged (%d, %d)\n", Nmodes, nconv);
    }

    /* Prepare headers in output files */
    if(me == 0) {
        opsec_info("Writing eigenmodes and eigenvalues to file...\n");

        Config opts = cfg_new();
        cfg_set(opts, "cellfile", cellfile);
        cfg_set_int(opts, "Ncells", Ncells);
        cfg_set_int(opts, "Nmodes", Nmodes);
        cfg_set_int(opts, "nconv", nconv);

        /* Write evalfile header */
        fprintf(feval, "# Signal matrix eigenvalues\n");
        abn_write_header(feval, nconv, REAL_FMT, opts);

        /* Write modefile header */
        fprintf(fmode, "# Signal matrix eigenmodes\n");
        fprintf(fmode, "# Each mode is stored consecutively, as below:\n");
        fprintf(fmode, "#   v_1,1   v_1,2   ... v_1,n\n");
        fprintf(fmode, "#   ...\n");
        fprintf(fmode, "#   v_m,1   v_m,2   ... v_m,n\n");
        abn_write_header(fmode, nconv*Ncells, REAL_FMT, opts);

        cfg_destroy(opts);
    }

    /* Create Ncells-by-1 vector r on root process, representing a single eigenvector */
    const slp::Context* pcontext = solver->GetContext();
    slp::Descriptor* rdesc = pcontext->new_descriptor(Ncells, 1, Ncells, 1);
    real* rvalues = (real*) opsec_malloc(rdesc->local_size() * sizeof(real));
    slp::Matrix<real> r(rdesc, &rvalues[0]);
    opsec_assert(rdesc->num_local_rows() == (me == 0 ? Ncells : 0));

    /* Gather eigenvectors to root process and write to file */
    for(int j = 0; j < nconv; j++) {
        real lambda = solver->GetEigenvalue(j);
        slp::Matrix<real> x = solver->GetEigenvector(j);

        slp::redistribute(Ncells, 1, x, 0, 0, r, 0, 0);

        if(me == 0) {
            feval.write((char*) &lambda, sizeof(real));
            fmode.write((char*) rvalues, Ncells*sizeof(real));
        }
    }

    /* Clean up nicely */
    if(me == 0) {
        feval.close();
        fmode.close();
    }
    free(rvalues);
    delete rdesc;
    delete solver;
    delete model;
    delete survey;
    free(cells);
    cfg_destroy(cfg);
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
