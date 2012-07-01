/* basis
 *
 * Break up the survey into a collection of basis cells in redshift space.
 * Use the selection function (provided by a Survey object) to compute the
 * expected number of galaxies in each cell in the absence of clustering.
 *
 * This program is designed to be run on a single node, using OpenMP for
 * shared-memory parallelism. */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "Cell.h"
#include "SelectionFunc.h"
#include "Survey.h"
#include "abn.h"
#include "cfg.h"
#include "cubature.h"
#include "opsec.h"
#include "rng.h"

template<class T> static inline T cb(T x) { return x*x*x; }

const char* usage =
    "Usage: %s [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h                Print this usage information\n"
    "  -c FILE           Read configuration options from FILE\n"
    "Configuration options:\n"
    "  coordsys=TYPE Coordinate system; either spherical or cartesian (required)\n"
    "  cellfile=FILE Output file for basis cells (required)\n"
    "  survey=NAME   Survey configuration (required)\n"
    "  Nx,Ny,Nz,XMin,XMax,YMin,YMax,ZMin,ZMax\n"
    "                Required if coordsys == 'cartesian'\n"
    "  Nr,Nmu,Nphi,RMin,RMax,MuMin,MuMax,PhiMin,PhiMax\n"
    "                Required if coordsys == 'spherical'\n";

struct CubatureData {
    SelectionFunc& nbar;
};

static void cartesian_integrand(unsigned ndim, unsigned npt, const double* uu,
                                void* data, unsigned fdim, double* fval)
{
    SelectionFunc& nbar = ((CubatureData*) data)->nbar;
    #pragma omp parallel for
    for(unsigned i = 0; i < npt; i++) {
        const double* u = &uu[3*i];
        double x = u[0];
        double y = u[1];
        double z = u[2];
        fval[i] = nbar(x, y, z);
    }
}

static void spherical_integrand(unsigned ndim, unsigned npt, const double* uu,
                                void* data, unsigned fdim, double* fval)
{
    SelectionFunc& nbar = ((CubatureData*) data)->nbar;
    #pragma omp parallel for
    for(unsigned i = 0; i < npt; i++) {
        const double* u = &uu[3*i];
        double r = cbrt(u[0]);
        double mu =     u[1];
        double phi =    u[2];
        fval[i] = nbar(r, mu, phi);
    }
}

/* Compute the effective volume and the expected number of galaxies within the
 * cell.  Return true if the cell should be kept, false if it should be
 * discarded. */
bool FinalizeCell(int coordsys, Cell& c, SelectionFunc& nbar, double epsrel = 1e-5, double epsabs = 1e-10) {
    /* Compute average of nbar within the cell */
    CubatureData data = { nbar };
    unsigned maxeval = 5000000;
    double N, err;
    if(coordsys == CoordSysCartesian) {
        double umin[3] = { c.xmin, c.ymin, c.zmin };
        double umax[3] = { c.xmax, c.ymax, c.zmax };
        adapt_integrate_v(1, cartesian_integrand, (void*) &data,
                          3, umin, umax,
                          maxeval, epsrel, epsabs,
                          &N, &err);
    }
    else if(coordsys == CoordSysSpherical) {
        double umin[3] = { cb(c.rmin), c.mumin, c.phimin };
        double umax[3] = { cb(c.rmax), c.mumax, c.phimax };
        adapt_integrate_v(1, spherical_integrand, (void*) &data,
                          3, umin, umax,
                          maxeval, epsrel, epsabs,
                          &N, &err);
    }

    /* Compute total volume of cell */
    double V = 0;
    if(coordsys == CoordSysCartesian) {
        V = (c.xmax - c.xmin) * (c.ymax - c.ymin) * (c.zmax - c.zmin);
    }
    else if(coordsys == CoordSysSpherical) {
        V = (cb(c.rmax) - cb(c.rmin))/3 * (c.mumax - c.mumin) * (c.phimax - c.phimin);
    }

    /* Keep cell only if nbar is not identically zero within the cell */
    if(N > 0) {
        c.Veff = V;
        c.Nbar = N;
        return true;
    }
    else
        return false;
}


int main(int argc, char* argv[]) {
    /* Initialize MPI, if available */
    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* This is a simple serial job, so only run on root process */
    if(me != 0) {
        MPI_Finalize();
        return 0;
    }

    /* Parse command line and configuration options */
    Config cfg = opsec_init(argc, argv, usage);

    /* Open output file first, to make sure we will be able to write */
    if(!cfg_has_key(cfg, "cellfile")) {
        opsec_error("basis: missing config option 'cellfile'\n");
        opsec_exit(1);
    }
    const char* cellfile = cfg_get(cfg, "cellfile");
    FILE* fcells = fopen(cellfile, "w");
    if(!fcells) {
        opsec_error("basis: cannot write to '%s'\n", cellfile);
        opsec_exit(1);
    }

    int coordsys = cfg_get_enum(cfg, "coordsys", "cartesian", CoordSysCartesian,
                                                 "spherical", CoordSysSpherical,
                                                 "", -1);
    if(coordsys == -1) {
        opsec_error("basis: invalid or missing config option: coordsys = %s\n", cfg_get(cfg, "coordsys"));
        opsec_exit(1);
    }

    /* Initialize survey */
    Survey* survey = InitializeSurvey(cfg);
    if(survey == NULL)
        opsec_exit(1);
    SelectionFunc nbar = survey->GetSelectionFunction();
    int N1 = survey->N1;
    int N2 = survey->N2;
    int N3 = survey->N3;

    vector<Cell> cells;
    int Ncells = 0;

    /* Populate list of cells */
    for(int d1 = 0; d1 < N1; d1++) {
        for(int d2 = 0; d2 < N2; d2++) {
            for(int d3 = 0; d3 < N3; d3++) {
                Cell c = survey->CreateEmptyCell(d1, d2, d3);

                if(FinalizeCell(coordsys, c, nbar)) {
                    c.a = Ncells++;
                    cells.push_back(c);
                }
            }
        }
    }

    opsec_info("Got %d non-empty cells out of %d possible.\n", Ncells, N1*N2*N3);

    Config opts = survey->GetConfigurationOptions();
    cfg_set_int(opts, "Ncells", Ncells);

    opsec_info("Writing cells to '%s'\n", cellfile);
    abn_write(fcells, &cells[0], Ncells, CELL_FMT_STRING, opts);

    fclose(fcells);
    cfg_destroy(opts);

    MPI_Finalize();
    return 0;
}
