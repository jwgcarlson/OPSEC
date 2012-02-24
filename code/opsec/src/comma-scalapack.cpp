/* comma-anasazi
 *
 * Compute the derivatives of the covariance matrix with respect to the model
 * parameters.  First compute this matrix in the cell basis, then use the KL
 * mode coefficients to project onto the KL basis.
 */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <string>
#include <vector>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "Cell.h"
#include "Model.h"
#include "SplitFile.h"
#include "Survey.h"
#include "XiFunc.h"
#include "abn.h"
#include "cfg.h"
#include "cscalapack.h"
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

/* BLACS context */
struct Context {
    int ictxt;      // BLACS context ID
    int nprow;      // number of rows in process grid
    int npcol;      // number of columns in process grid
    int myrow;      // my row
    int mycol;      // my column
};

Context create_context(int numrows, int numcols) {
    Context context;
    Cblacs_get(0, 0, &context.ictxt);
    Cblacs_gridinit(&context.ictxt, "Row", numrows, numcols);
    Cblacs_gridinfo(context.ictxt, &context.nprow, &context.npcol, &context.myrow, &context.mycol);
    if(context.myrow >= numrows || context.mycol >= numcols)
        context.myrow = context.mycol = -1;
    return context;
}

/* BLACS array descriptor */
struct ArrayDesc {
    int dtype;  // = 1 for dense matrices
    int ictxt;  // BLACS context handle
    int m;      // number of rows
    int n;      // number of columns
    int mb;     // row blocking factor
    int nb;     // column blocking factor
    int rsrc;   // = 0
    int csrc;   // = 0
    int lld;    // leading dimension of local array
};

static inline int mynumroc(int n, int nb, int iproc, int isrcproc, int nprocs) {
    return numroc_(&n, &nb, &iproc, &isrcproc, &nprocs);
}

static ArrayDesc create_descriptor(const Context& context,
                                   int m, int n, int mb, int nb,
                                   int rsrc = 0, int csrc = 0, int lld = 0)
{
    ArrayDesc desc;
    desc.dtype = 1;
    desc.ictxt = context.ictxt;
    desc.m = m;
    desc.n = n;
    desc.mb = mb;
    desc.nb = nb;
    desc.rsrc = rsrc;
    desc.csrc = csrc;
    desc.lld = std::max(lld, mynumroc(m, mb, context.myrow, rsrc, context.nprow));
    return desc;
}

static inline int myindxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs) {
    iloc += 1;
    return indxl2g_(&iloc, &nb, &iproc, &isrcproc, &nprocs) - 1;
}

static inline void mypdgemv(const char* transa, int m, int n,
                            double alpha,
                            const double* a, int ia, int ja, const ArrayDesc& desca,
                            const double* x, int ix, int jx, const ArrayDesc& descx, int incx,
                            double beta,
                            double* y, int iy, int jy, const ArrayDesc& descy, int incy)
{
    ia++; ja++; ix++; jx++; iy++; jy++;         // convert to Fortran indexing
    pdgemv_((char*) transa, &m, &n,
            &alpha,
            (double*) a, &ia, &ja, (int*) &desca,
            (double*) x, &ix, &jx, (int*) &descx, &incx,
            &beta,
            y, &iy, &jy, (int*) &descy, &incy);
}

static inline void mypdgemm(const char* transa, const char* transb, int m, int n, int k,
                            double alpha,
                            const double* a, int ia, int ja, const ArrayDesc& desca,
                            const double* b, int ib, int jb, const ArrayDesc& descb,
                            double beta,
                            double* c, int ic, int jc, const ArrayDesc& descc)
{
    ia++; ja++; ib++; jb++; ic++; jc++;         // convert to Fortran indexing
    pdgemm_((char*) transa, (char*) transb, &m, &n, &k,
            &alpha,
            (double*) a, &ia, &ja, (int*) &desca,
            (double*) b, &ib, &jb, (int*) &descb,
            &beta,
            (double*) c, &ic, &jc, (int*) &descc);
}

extern "C" void Cpdgemr2d(int m, int n, double* a, int ia, int ja, int* desca, double* b, int ib, int jb, int* descb, int ictxt);
static inline void mypdgemr2d(int m, int n,
                              const double* a, int ia, int ja, const ArrayDesc& desca,
                              double* b, int ib, int jb, const ArrayDesc& descb,
                              const Context& gcontext)
{
    Cpdgemr2d(m, n,
              (double*) a, ia + 1, ja + 1, (int*) &desca, 
              b, ib + 1, jb + 1, (int*) &descb, 
              gcontext.ictxt);
}

int main(int argc, char* argv[]) {
    int Ncells, Nmodes, Nparams;
    const char *cellfile, *modefile, *covfile;

#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    /* Set up global BLACS context encompassing all processes */
    int nprocs, me;
    Cblacs_pinfo(&me, &nprocs);
    Context gcontext = create_context(nprocs, 1);

    Config cfg = opsec_init(argc, argv, usage);

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

    /* Decide how best to partition available processes into 2-D process grid */
    int nprow = nprocs, npcol = 1;
    if(cfg_has_keys(cfg, "nprow,npcol", ",")) {
        /* Use explicit values */
        nprow = cfg_get_int(cfg, "nprow");
        npcol = cfg_get_int(cfg, "npcol");
        if(nprow*npcol != nprocs) {
            fprintf(stderr, "comma: requested %d-by-%d process grid on %d processes\n", nprow, npcol, nprocs);
            opsec_exit(1);
        }
    }
    else {
        /* Make the process grid square-ish if possible */
        while(nprow > npcol && (nprow % 2) == 0) {
            nprow /= 2;
            npcol *= 2;
        }
    }
    if(me == 0)
        printf("Using %d-by-%d process grid\n", nprow, npcol);

    /* Set up a new BLACS context for this 2-D process grid */
    Context pcontext = create_context(nprow, npcol);

    /* Set up a third BLACS context (on the root process only) for data I/O */
    Context rcontext = create_context(1, 1);
    if(me == 0)
        assert(rcontext.myrow == 0 && rcontext.mycol == 0);

    /* Determine coordinate system */
    int coordsys;
    const char* coordsysstr = cfg_get(cfg, "coordsys");
    if(strcmp(coordsysstr, "spherical") == 0)
        coordsys = CoordSysSpherical;
    else if(strcmp(coordsysstr, "cartesian") == 0)
        coordsys = CoordSysCartesian;
    else {
        fprintf(stderr, "comma: missing or invalid config option: coordsys = %s\n", coordsysstr);
        opsec_exit(1);
    }

    int N1, N2, N3;
    if(coordsys == CoordSysSpherical) {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi", ",")) {
            if(me == 0)
                fprintf(stderr, "comma: need Nr,Nmu,Nphi for spherical coordinates\n");
            opsec_exit(1);
        }
        N1 = cfg_get_int(cfg, "Nr");
        N2 = cfg_get_int(cfg, "Nmu");
        N3 = cfg_get_int(cfg, "Nphi");
    }
    else if(coordsys == CoordSysCartesian) {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz", ",")) {
            if(me == 0)
                fprintf(stderr, "comma: need Nx,Ny,Nz for cartesian coordinates\n");
            opsec_exit(1);
        }
        N1 = cfg_get_int(cfg, "Nx");
        N2 = cfg_get_int(cfg, "Ny");
        N3 = cfg_get_int(cfg, "Nz");
    }

    /* Load model */
    Model* model = InitializeModel(cfg);
    if(!model)
        opsec_exit(1);
    Nparams = model->NumParams();

    /* Load survey */
    Survey* survey = InitializeSurvey(cfg);
    if(!survey)
        opsec_exit(1);

    /* Read cells from file */
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);

    /* Define array descriptor for the mode coefficient matrix B */
    int m_B = Ncells;
    int n_B = Nmodes;
    int mb_B = std::min(32, m_B/nprocs);
    int nb_B = std::min(32, n_B/nprocs);
    int mloc_B = mynumroc(m_B, mb_B, pcontext.myrow, 0, pcontext.nprow);
    int nloc_B = mynumroc(n_B, nb_B, pcontext.mycol, 0, pcontext.npcol);
    ArrayDesc descB = create_descriptor(pcontext, m_B, n_B, mb_B, nb_B, 0, 0, mloc_B);
    real* Bloc = (real*) opsec_malloc(mloc_B*nloc_B*sizeof(real));

    /* Define array descriptor for a single mode */
    ArrayDesc descmode = create_descriptor(rcontext, m_B, 1, m_B, 1, 0, 0, m_B);
    real* mode = NULL;
    if(me == 0)
        mode = (real*) opsec_malloc(m_B*sizeof(real));

    /* Read modes file header */
    SplitFile fmodes;
    if(me == 0) {
        fmodes.open(modefile, "r");
        if(!fmodes.isopen())
            opsec_abort(1);

        size_t n, size;
        char endian, fmt[ABN_MAX_FORMAT_LENGTH];
        Config opts = cfg_new();
        abn_read_header(fmodes, &n, &size, &endian, fmt, opts);
        if(cfg_has_key(opts, "Nmodes"))
            assert(cfg_get_int(opts, "Nmodes") == Nmodes);
        if(cfg_has_key(opts, "Ncells"))
            assert(cfg_get_int(opts, "Ncells") == Ncells);
        cfg_destroy(opts);
    }

    /* Read modes one at a time, redistributing from root process to full process grid */
    for(int j = 0; j < Nmodes; j++) {
        if(me == 0) {
            ssize_t nread = fmodes.read((char*) mode, Ncells*sizeof(real));
            if(nread != Ncells*sizeof(real)) {
                fprintf(stderr, "comma: error reading modes (%s,%d,%zd,%d)\n", fmodes.get_filename().c_str(), j, nread, Ncells);
                opsec_abort(1);
            }
        }

        mypdgemr2d(m_B, 1, mode, 0, 0, descmode, Bloc, 0, j, descB, gcontext);
    }

    if(me == 0) {
        fmodes.close();
        free(mode);
    }

    /* Define array descriptor for the signal matrix S in the cell basis */
    int m_S = Ncells;
    int n_S = Ncells;
    int mb_S = std::min(32, m_S/nprocs);
    int nb_S = std::min(32, n_S/nprocs);
    int mloc_S = mynumroc(m_S, mb_S, pcontext.myrow, 0, pcontext.nprow);
    int nloc_S = mynumroc(n_S, nb_S, pcontext.mycol, 0, pcontext.npcol);
    ArrayDesc descS = create_descriptor(pcontext, m_S, n_S, mb_S, nb_S, 0, 0, mloc_S);
    real* Sloc = (real*) opsec_malloc(mloc_S*nloc_S*sizeof(real));

    /* Define array descriptor for the product matrix C = S . B */
    int m_C = Ncells;
    int n_C = Nmodes;
    int mb_C = std::min(32, m_C/nprocs);
    int nb_C = std::min(32, n_C/nprocs);
    int mloc_C = mynumroc(m_C, mb_C, pcontext.myrow, 0, pcontext.nprow);
    int nloc_C = mynumroc(n_C, nb_C, pcontext.mycol, 0, pcontext.npcol);
    ArrayDesc descC = create_descriptor(pcontext, m_C, n_C, mb_C, nb_C, 0, 0, mloc_C);
    real* Cloc = (real*) opsec_malloc(mloc_C*nloc_C*sizeof(real));

    /* Define array descriptor for the projected signal matrix P = B^T . S . B = B^T . C */
    int m_P = Nmodes;
    int n_P = Nmodes;
    int mb_P = std::min(32, m_P/nprocs);
    int nb_P = std::min(32, n_P/nprocs);
    int mloc_P = mynumroc(m_P, mb_P, pcontext.myrow, 0, pcontext.nprow);
    int nloc_P = mynumroc(n_P, nb_P, pcontext.mycol, 0, pcontext.npcol);
    ArrayDesc descP = create_descriptor(pcontext, m_P, n_P, mb_P, nb_P, 0, 0, mloc_P);
    /* P will reuse the same memory as S, since they never need to be used at the same time */
    real* Ploc = Sloc;

    /* Define array descriptor for a single column of P (on the root process) */
    ArrayDesc desccolP = create_descriptor(rcontext, m_P, 1, m_P, 1, 0, 0, m_P);
    real* colP = NULL;
    if(me == 0)
        colP = (real*) opsec_malloc(m_P*sizeof(real));

    /* Prepare output file */
    SplitFile fout;
    if(me == 0) {
        fout.open(covfile, "w");
        if(!fout.isopen()) {
            fprintf(stderr, "comma: could not open '%s' for writing\n", covfile);
            opsec_exit(1);
        }

        fprintf(fout, "# Derivatives of the covariance matrix in KL basis, C_{ij,m}.\n");
        fprintf(fout, "# Binary data consists of %d consecutive symmetric %dx%d matrices\n", Nparams, Nmodes, Nmodes);
        fprintf(fout, "# in packed storage.  Each matrix consists of %d elements.\n", Nmodes*(Nmodes+1)/2);

        Config opts = cfg_new();
        cfg_set_int(opts, "Nmodes", Nmodes);
        cfg_set_int(opts, "Nparams", Nparams);
        abn_write_header(fout, Nparams*Nmodes*(Nmodes+1)/2, "d", opts);
        cfg_destroy(opts);
    }

    for(int param = 0; param < Nparams; param++) {
        if(me == 0) {
            printf("Parameter %d/%d\n", param+1, Nparams);
            fflush(stdout);
        }

        XiFunc xi = model->GetXiDeriv(param);

        if(me == 0) {
            printf("  Computing matrix elements of C,%d in cell basis\n", param+1);
            fflush(stdout);
        }

        /* Compute local signal matrix elements in cell basis */
        std::vector<int> rows(mloc_S), cols(nloc_S);
        for(int aloc = 0; aloc < mloc_S; aloc++)
            rows[aloc] = myindxl2g(aloc, descS.mb, pcontext.myrow, descS.rsrc, pcontext.nprow);
        for(int bloc = 0; bloc < nloc_S; bloc++)
            cols[bloc] = myindxl2g(bloc, descS.nb, pcontext.mycol, descS.csrc, pcontext.npcol);
        ComputeSignalMatrixBlock(Ncells, mloc_S, &rows[0], nloc_S, &cols[0],
                                 Sloc, descS.lld, xi, survey, coordsys, cells, N1, N2, N3);

        if(me == 0) {
            printf("  Projecting onto mode basis\n");
            fflush(stdout);
        }

        /* Compute C = S . B */
        mypdgemm("N", "N", Ncells, Nmodes, Ncells, 1.0, Sloc, 0, 0, descS, Bloc, 0, 0, descB, 0.0, Cloc, 0, 0, descC);

        /* Compute P = B^T . S */
        mypdgemm("T", "N", Nmodes, Nmodes, Ncells, 1.0, Bloc, 0, 0, descB, Cloc, 0, 0, descC, 0.0, Ploc, 0, 0, descP);

        if(me == 0) {
            printf("  Gathering results on root process and writing to file\n");
            fflush(stdout);
        }

        /* Gather P to root process, column by column, and write to file */
        for(int j = 0; j < Nmodes; j++) {
            mypdgemr2d(m_P, 1, Ploc, 0, j, descP, colP, 0, 0, desccolP, gcontext);
            if(me == 0) {
                fout.write((char*) colP, (j+1)*sizeof(real));
                for(int i = 0; i <= j; i++)
                    printf("D[%d](%d,%d) = %g\n", param, i, j, colP[i]);
            }
        }
    }
    if(me == 0)
        fout.close();

    /* Clean up nicely */
    free(cells);
    free(Bloc);
    free(Cloc);
    free(Sloc);
//    Do not free(Ploc) since Ploc = Sloc
    free(colP);
    delete model;
    Cblacs_gridexit(rcontext.ictxt);
    Cblacs_gridexit(pcontext.ictxt);
    Cblacs_gridexit(gcontext.ictxt);
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
