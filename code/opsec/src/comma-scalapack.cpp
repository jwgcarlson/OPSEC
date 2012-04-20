/* comma
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
#include "slp.h"

const char* usage =
    "Usage: comma [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h             Display this help\n"
    "  -c FILE        Read additional configuration options from FILE\n"
    "Configuration options:\n"
    "  cellfile=FILE  Read list of cells from FILE (required)\n"
    "  modefile=FILE  Read KL modes from FILE (required)\n"
    "  covfile=FILE   Write covariance matrices to FILE (required)\n"
    "  Nmodes=NUM     Number of KL modes (required)\n"
    "  model=NAME     Use specified model (required)\n"
    "  coordsys=TYPE  Coordinate system; either spherical or cartesian (required)\n"
    "  Nr,Nmu,Nphi    Spherical grid (required if coordsys = spherical)\n"
    "  Nx,Ny,Nz       Cartesian grid (required if coordsys = cartesian)\n";

/* Read mode matrix B from file */
void ReadModes(const slp::Context& pcontext, const char* modefile, int Ncells, int Nmodes, slp::Matrix<real>& B) {
    int me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* Create Ncells-by-1 vector on the root process */
    slp::Descriptor* zdesc = pcontext.new_descriptor(Ncells, 1, Ncells, 1);
    real* zvalues = (real*) opsec_malloc(zdesc->local_size() * sizeof(real));
    slp::Matrix<real> z(zdesc, zvalues);

    /* Read modes file header */
    SplitFile fmodes;
    if(me == 0) {
        fmodes.open(modefile, "r");
        if(!fmodes.isopen()) {
            opsec_error("comma: cannot open modes file %s\n", modefile);
            opsec_abort(1);
        }

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
            size_t nread = fmodes.read((char*) zvalues, Ncells*sizeof(real));
            if(nread != Ncells*sizeof(real)) {
                opsec_error("comma: error reading modes (%s,%d,%zd,%d)\n", fmodes.get_filename().c_str(), j, nread, Ncells);
                opsec_abort(1);
            }
        }

        slp::redistribute(Ncells, 1, z, 0, 0, B, 0, j);
    }

    delete zdesc;
    free(zvalues);
}

int main(int argc, char* argv[]) {
    /* Initialize MPI */
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* Parse command line and configuration options */
    Config cfg = opsec_init(argc, argv, usage);

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "coordsys,cellfile,modefile,covfile,Nmodes", ",")) {
        if(me == 0) opsec_error("comma: missing configuration options\n");
        if(me == 0) fputs(usage, stderr);
        opsec_exit(1);
    }
    const char* cellfile = cfg_get(cfg, "cellfile");
    const char* modefile = cfg_get(cfg, "modefile");
    const char* covfile = cfg_get(cfg, "covfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");

    /* Decide how best to partition available processes into 2-D process grid */
    int nprow, npcol;
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
        nprow = nprocs;
        npcol = 1;
        while(nprow > npcol && (nprow % 2) == 0) {
            nprow /= 2;
            npcol *= 2;
        }
    }
    if(me == 0)
        opsec_info("Using %d-by-%d process grid\n", nprow, npcol);

    /* Set up a new BLACS context for this 2-D process grid */
    slp::Context pcontext(nprow, npcol);

    /* Determine coordinate system */
    int coordsys = cfg_get_enum(cfg, "coordsys", "cartesian", CoordSysCartesian,
                                                 "spherical", CoordSysSpherical,
                                                 "", -1);
    if(coordsys == -1) {
        if(me == 0)
            opsec_error("comma: missing or invalid config option: coordsys = %s\n", cfg_get(cfg, "coordsys"));
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
    int Nparams = model->NumParams();

    /* Load survey */
    Survey* survey = InitializeSurvey(cfg);
    if(!survey)
        opsec_exit(1);

    /* Read cells from file */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL)
        opsec_exit(1);
    assert(Ncells >= Nmodes);

    /* Choose reasonable blocking factors */
    int Ncells_block, Nmodes_block;
    if(nprocs == 1) {
        Ncells_block = Ncells;
        Nmodes_block = Nmodes;
    }
    else {
        Ncells_block = std::min(Ncells/nprow, 64);
        Nmodes_block = std::min(Nmodes/npcol, 64);
    }

    /* Create Ncells-by-Nmodes mode coefficient matrix B */
    slp::Descriptor* Bdesc = pcontext.new_descriptor(Ncells, Nmodes, Ncells_block, Nmodes_block);
    real* Bvalues = (real*) opsec_malloc(Bdesc->local_size() * sizeof(real));
    slp::Matrix<real> B(Bdesc, Bvalues);

    /* Read in B */
    ReadModes(pcontext, modefile, Ncells, Nmodes, B);

    /* Create Ncells-by-Ncells signal matrix S */
    slp::Descriptor* Sdesc = pcontext.new_descriptor(Ncells, Ncells, Ncells_block, Ncells_block);
    real* Svalues = (real*) opsec_malloc(Sdesc->local_size() * sizeof(real));
    slp::Matrix<real> S(Sdesc, Svalues);

    /* Create Ncells-by-Nmodes product matrix C = S*B */
    slp::Descriptor* Cdesc = pcontext.new_descriptor(Ncells, Nmodes, Ncells_block, Nmodes_block);
    real* Cvalues = (real*) opsec_malloc(Cdesc->local_size() * sizeof(real));
    slp::Matrix<real> C(Cdesc, Cvalues);

    /* Create Nmodes-by-Nmodes projected signal matrix P = B^T * S * B = B^T * C */
    slp::Descriptor* Pdesc = pcontext.new_descriptor(Nmodes, Nmodes, Nmodes_block, Nmodes_block);
    real* Pvalues = Svalues;    // P can use same memory as S, since there is no overlap
    slp::Matrix<real> P(Pdesc, Pvalues);

    /* Create Nmodes-by-1 vector r on root process, representing a single column of P */
    slp::Descriptor* rdesc = pcontext.new_descriptor(Nmodes, 1, Nmodes, 1);
    real* rvalues = (real*) opsec_malloc(rdesc->local_size() * sizeof(real));
    slp::Matrix<real> r(rdesc, rvalues);

    /* Prepare output file */
    SplitFile fout;
    if(me == 0) {
        fout.open(covfile, "w");
        if(!fout.isopen()) {
            opsec_error("comma: could not open '%s' for writing\n", covfile);
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
        if(me == 0)
            opsec_info("Parameter %d/%d\n", param+1, Nparams);

        XiFunc xi = model->GetXiDeriv(param);

        if(me == 0)
            opsec_info("  Computing matrix elements of C,%d in cell basis\n", param+1);

        /* Compute local signal matrix elements in cell basis */
        std::vector<int> rows(Sdesc->mloc), cols(Sdesc->nloc);
        for(int aloc = 0; aloc < Sdesc->mloc; aloc++)
            rows[aloc] = Sdesc->row_l2g(aloc);
        for(int bloc = 0; bloc < Sdesc->nloc; bloc++)
            cols[bloc] = Sdesc->col_l2g(bloc);
        ComputeSignalMatrixBlock(Ncells, rows, cols, Svalues, Sdesc->lld,
                                 xi, survey, cells, N1, N2, N3);

        if(me == 0)
            opsec_info("  Projecting onto mode basis\n");

        /* Compute C = S * B */
        slp::multiply(S, B, C);

        /* Compute P = B^T * C */
        slp::multiply(B, C, P, 'T', 'N');

        if(me == 0)
            opsec_info("  Gathering results on root process and writing to file\n");

        /* Gather P to root process, column by column, and write to file */
        for(int j = 0; j < Nmodes; j++) {
            slp::redistribute(Nmodes, 1, P, 0, j, r, 0, 0);
            if(me == 0) {
                fout.write((char*) rvalues, (j+1)*sizeof(real));
//                for(int i = 0; i <= j; i++)
//                    printf("D[%d](%d,%d) = %g\n", param, i, j, rvalues[i]);
            }
        }
    }
    if(me == 0)
        fout.close();

    /* Clean up nicely */
    free(cells);
    free(Bvalues);
    free(Cvalues);
    free(Svalues);
//    Do not free(Pvalues) since Pvalues = Svalues
    free(rvalues);
    delete Bdesc;
    delete Cdesc;
    delete Sdesc;
    delete Pdesc;
    delete rdesc;
    pcontext.exit();
    delete model;
#ifdef OPSEC_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
