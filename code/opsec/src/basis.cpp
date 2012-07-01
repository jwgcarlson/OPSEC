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

    /* Initialize selection function */
    Survey* survey = InitializeSurvey(cfg);
    if(survey == NULL)
        opsec_exit(1);
    SelectionFunc nbar = survey->GetSelectionFunction();

    vector<Cell> cells;
    int Ncells = 0;

    int N1, N2, N3;
    double Min1, Max1, Min2, Max2, Min3, Max3;
    Config opts = cfg_new();
    if(coordsys == CoordSysSpherical) {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi,RMin,RMax,MuMin,MuMax,PhiMin,PhiMax", ",")) {
            opsec_error("basis: must provide config options N{r,mu,phi} and {R,Mu,Phi}{Min,Max}\n");
            opsec_exit(1);
        }

        N1 = cfg_get_int(cfg, "Nr");
        N2 = cfg_get_int(cfg, "Nmu");
        N3 = cfg_get_int(cfg, "Nphi");
        Min1 = cfg_get_double(cfg, "RMin");
        Max1 = cfg_get_double(cfg, "RMax");
        Min2 = cfg_get_double(cfg, "MuMin");
        Max2 = cfg_get_double(cfg, "MuMax");
        Min3 = cfg_get_double(cfg, "PhiMin");
        Max3 = cfg_get_double(cfg, "PhiMax");

        fprintf(fcells, "# struct Cell {\n");
        fprintf(fcells, "#     int a;\n");
        fprintf(fcells, "#     int G;\n");
        fprintf(fcells, "#     double rmin, rmax;\n");
        fprintf(fcells, "#     double mumin, mumax;\n");
        fprintf(fcells, "#     double phimin, phimax;\n");
        fprintf(fcells, "#     double Veff;\n");
        fprintf(fcells, "#     double Nbar;\n");
        fprintf(fcells, "# };\n");

        cfg_set(opts, "coordsys", "spherical");
        cfg_set_int(opts, "Nr", N1);
        cfg_set_int(opts, "Nmu", N2);
        cfg_set_int(opts, "Nphi", N3);
        cfg_set_double(opts, "RMin", Min1);
        cfg_set_double(opts, "RMax", Max1);
        cfg_set_double(opts, "MuMin", Min2);
        cfg_set_double(opts, "MuMax", Max2);
        cfg_set_double(opts, "PhiMin", Min3);
        cfg_set_double(opts, "PhiMax", Max3);
    }
    else if(coordsys == CoordSysCartesian) {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz,XMin,XMax,YMin,YMax,ZMin,ZMax", ",")) {
            opsec_error("basis: must provide config options N{x,y,z} and {X,Y,Phi}{Z,Max}\n");
            opsec_exit(1);
        }

        N1 = cfg_get_int(cfg, "Nx");
        N2 = cfg_get_int(cfg, "Ny");
        N3 = cfg_get_int(cfg, "Nz");
        Min1 = cfg_get_double(cfg, "XMin");
        Max1 = cfg_get_double(cfg, "XMax");
        Min2 = cfg_get_double(cfg, "YMin");
        Max2 = cfg_get_double(cfg, "YMax");
        Min3 = cfg_get_double(cfg, "ZMin");
        Max3 = cfg_get_double(cfg, "ZMax");

        fprintf(fcells, "# struct Cell {\n");
        fprintf(fcells, "#     int a;\n");
        fprintf(fcells, "#     int G;\n");
        fprintf(fcells, "#     double xmin, xmax;\n");
        fprintf(fcells, "#     double ymin, ymax;\n");
        fprintf(fcells, "#     double zmin, zmax;\n");
        fprintf(fcells, "#     double Veff;\n");
        fprintf(fcells, "#     double Nbar;\n");
        fprintf(fcells, "# };\n");


        cfg_set(opts, "coordsys", "cartesian");
        cfg_set_int(opts, "Nx", N1);
        cfg_set_int(opts, "Ny", N2);
        cfg_set_int(opts, "Nz", N3);
        cfg_set_double(opts, "XMin", Min1);
        cfg_set_double(opts, "XMax", Max1);
        cfg_set_double(opts, "YMin", Min2);
        cfg_set_double(opts, "YMax", Max2);
        cfg_set_double(opts, "ZMin", Min3);
        cfg_set_double(opts, "ZMax", Max3);
    }

    /* Populate list of cells */
    Cell c;
    #pragma omp parallel for private(c)
    for(int d = 0; d < N1; d++) {
        c.min1 = Min1 +     d*(Max1 - Min1)/N1;
        c.max1 = Min1 + (d+1)*(Max1 - Min1)/N1;
        for(int e = 0; e < N2; e++) {
            c.min2 = Min2 +     e*(Max2 - Min2)/N2;
            c.max2 = Min2 + (e+1)*(Max2 - Min2)/N2;
            for(int f = 0; f < N3; f++) {
                c.min3 = Min3 +     f*(Max3 - Min3)/N3;
                c.max3 = Min3 + (f+1)*(Max3 - Min3)/N3;
                c.G = N3*(N2*d + e) + f;

                if(FinalizeCell(coordsys, c, nbar)) {
                    #pragma omp critical (cells_update)
                    {
                        c.a = Ncells++;
                        cells.push_back(c);
                    }
                }
            }
        }
    }

    opsec_info("Got %d non-empty cells out of %d possible.\n", Ncells, N1*N2*N3);
    cfg_set_int(opts, "Ncells", Ncells);

    opsec_info("Writing cells to '%s'.\n", cellfile);
    abn_write(fcells, &cells[0], Ncells, CELL_FMT_STRING, opts);

    fclose(fcells);
    cfg_destroy(opts);

    MPI_Finalize();
    return 0;
}
