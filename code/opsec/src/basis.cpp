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

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "Cell.h"
#include "Survey.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "rng.h"

template<typename T> static inline T cb(T x) { return x*x*x; }

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

#if 0
void exit_basis(int status) {
    exit(status);
}
#endif

/* Compute the effective volume of and the expected number of galaxies within
 * the cell.  Return true if the cell should be kept, false if it should be
 * discarded. */
bool FinalizeCellS(Cell& c, SelectionFunc& nbar, double epsrel = 1e-5, double epsabs = 1e-10) {
    const int lg2n_min = 6, lg2n_max = 30;

    Sobol sobol;
    rng_sobol_init(&sobol, 3);
    double p[3];

    double r, mu, phi;
    double rmincb = cb(c.rmin), rmaxcb = cb(c.rmax);

    int n = 0, numzero = 0;
    double f, fsum = 0.;
    double Q, oldQ = 0., dQ;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        while(n < (1 << lg2n)) {
            rng_sobol_get(&sobol, p);

            r = cbrt((1 - p[0]) * rmincb + p[0] * rmaxcb);
            mu = (1 - p[1]) * c.mumin + p[1] * c.mumax;
            phi = (1 - p[2]) * c.phimin + p[2] * c.phimax;
            f = nbar(r, mu, phi);
            if(f == 0.0)
                numzero++;
            else
                fsum += f;
            n++;
        }

        /* Declare convergence if the integral changes by only a small amount
         * after doubling the number of evaluation points. */
        Q = fsum/n;
        dQ = fabs(Q - oldQ);
        if(dQ < fmax(epsabs, fabs(Q)*epsrel))
            break;

        /* Otherwise continue with the next batch of evaluations. */
        oldQ = Q;
    }

    /* Total volume of cell */
    double V = (rmaxcb - rmincb)/3 * (c.mumax - c.mumin) * (c.phimax - c.phimin);

    /* Fraction of cell for which the number density vanishes */
    double fzero = double(numzero)/double(n);

    /* Keep cell only if nbar is nonzero over a significant fraction of the volume */
    if(fzero > 0.2)
        return false;
    else {
        c.Veff = (1 - fzero) * V;
        c.Nbar = Q * V;
        return true;
    }
}

/* Compute the effective volume of and the expected number of galaxies within
 * the cell.  Return true if the cell should be kept, false if it should be
 * discarded. */
bool FinalizeCellC(Cell& c, SelectionFunc& nbar, double epsrel = 1e-5, double epsabs = 1e-10) {
    const int lg2n_min = 6, lg2n_max = 30;

    Sobol sobol;
    rng_sobol_init(&sobol, 3);
    double p[3];

    double x, y, z;
    int n = 0, numzero = 0;
    double f, fsum = 0.;
    double Q, oldQ = 0., dQ;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        while(n < (1 << lg2n)) {
            rng_sobol_get(&sobol, p);

            x = (1 - p[0])*c.xmin + p[0]*c.xmax;
            y = (1 - p[1])*c.ymin + p[1]*c.ymax;
            z = (1 - p[2])*c.zmin + p[2]*c.zmax;
            f = nbar(x, y, z);
            if(f == 0.0)
                numzero++;
            else
                fsum += f;
            n++;
        }

        /* Declare convergence if the integral changes by only a small amount
         * after doubling the number of evaluation points. */
        Q = fsum/n;
        dQ = fabs(Q - oldQ);
        if(dQ < fmax(epsabs, fabs(Q)*epsrel))
            break;

        /* Otherwise continue with the next batch of evaluations. */
        oldQ = Q;
    }

    /* Total volume of cell */
    double V = (c.xmax - c.xmin) * (c.ymax - c.ymin) * (c.zmax - c.zmin);

    /* Fraction of cell for which the number density vanishes */
    double fzero = double(numzero)/double(n);

    /* Keep cell only if nbar is nonzero over a significant fraction of the volume */
    if(fzero > 0.2)
        return false;
    else {
        c.Veff = (1 - fzero) * V;
        c.Nbar = Q * V;
        return true;
    }
}


int main(int argc, char* argv[]) {
    Config cfg = opsec_init(argc, argv, usage);
#if 0
    Config cfg = cfg_new();

    /* Parse command line switches */
    int opt;
    const char* optstring = "hc:";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            fputs(usage, stdout);
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fputs(usage, stderr);
            return 1;
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);
#endif

    /* Debugging... */
//    printf("# Config options\n");
//    cfg_write(cfg, stdout);
//    printf("\n");

    /* Open output file first, to make sure we will be able to write */
    if(!cfg_has_key(cfg, "cellfile")) {
        fprintf(stderr, "basis: missing config option 'cellfile'\n");
//        return 1;
        opsec_exit(1);
    }
    const char* cellfile = cfg_get(cfg, "cellfile");
    FILE* fcells = fopen(cellfile, "w");
    if(!fcells) {
        fprintf(stderr, "basis: cannot write cells to file '%s'\n", cellfile);
//        return 1;
        opsec_exit(1);
    }

    string coordsys = cfg_get(cfg, "coordsys");
    if(coordsys != "spherical" && coordsys != "cartesian") {
        fprintf(stderr, "basis: invalid or missing config option 'coordsys'\n");
//        return 1;
        opsec_exit(1);
    }

    /* Initialize selection function */
    Survey* survey = InitializeSurvey(cfg);
    if(survey == NULL)
//        exit_basis(1);
        opsec_exit(1);
    SelectionFunc nbar = survey->GetSelectionFunction();

    vector<Cell> cells;
    int Ncells = 0;

    if(coordsys == "spherical") {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi,RMin,RMax,MuMin,MuMax,PhiMin,PhiMax", ",")) {
            fprintf(stderr, "basis: must provide config options N{r,mu,phi} and {R,Mu,Phi}{Min,Max}\n");
//            return 1;
            opsec_exit(1);
        }

        int Nr = cfg_get_int(cfg, "Nr");
        int Nmu = cfg_get_int(cfg, "Nmu");
        int Nphi = cfg_get_int(cfg, "Nphi");
        double RMin = cfg_get_double(cfg, "RMin");
        double RMax = cfg_get_double(cfg, "RMax");
        double MuMin = cfg_get_double(cfg, "MuMin");
        double MuMax = cfg_get_double(cfg, "MuMax");
        double PhiMin = cfg_get_double(cfg, "PhiMin");
        double PhiMax = cfg_get_double(cfg, "PhiMax");

        /* Populate list of cells */
        Cell c;
        int d, e, f;
        #pragma omp parallel for private(c,d,e,f)
        for(d = 0; d < Nr; d++) {
            c.rmin = RMin +     d*(RMax - RMin)/Nr;   // this choice of linear radial spacing is consistent with the assumptions in dot.cpp
            c.rmax = RMin + (d+1)*(RMax - RMin)/Nr;
            for(e = 0; e < Nmu; e++) {
                c.mumin = MuMin +     e*(MuMax - MuMin)/Nmu;
                c.mumax = MuMin + (e+1)*(MuMax - MuMin)/Nmu;
                for(f = 0; f < Nphi; f++) {
                    c.phimin = PhiMin +     f*(PhiMax - PhiMin)/Nphi;
                    c.phimax = PhiMin + (f+1)*(PhiMax - PhiMin)/Nphi;
                    c.G = f + Nphi*e + Nmu*Nphi*d;

                    if(FinalizeCellS(c, nbar, 1e-4, 1e-12)) {
                        #pragma omp critical (cells_update)
                        {
                            c.a = Ncells++;
                            cells.push_back(c);
                        }
                    }
                }
            }
        }

        printf("Got %d non-empty cells out of %d possible.\n", Ncells, Nr*Nmu*Nphi);
        printf("Writing list of cells to '%s'.\n", cellfile);

        fprintf(fcells, "# struct Cell {\n");
        fprintf(fcells, "#     int a;\n");
        fprintf(fcells, "#     int G;\n");
        fprintf(fcells, "#     double rmin, rmax;\n");
        fprintf(fcells, "#     double mumin, mumax;\n");
        fprintf(fcells, "#     double phimin, phimax;\n");
        fprintf(fcells, "#     double Veff;\n");
        fprintf(fcells, "#     double Nbar;\n");
        fprintf(fcells, "# };\n");

        Config opts = cfg_new();
        cfg_set(opts, "coordsys", "spherical");
        cfg_set_int(opts, "Ncells", Ncells);
        cfg_set_int(opts, "Nr", Nr);
        cfg_set_int(opts, "Nmu", Nmu);
        cfg_set_int(opts, "Nphi", Nphi);
        cfg_set_double(opts, "RMin", RMin);
        cfg_set_double(opts, "RMax", RMax);
        cfg_set_double(opts, "MuMin", MuMin);
        cfg_set_double(opts, "MuMax", MuMax);
        cfg_set_double(opts, "PhiMin", PhiMin);
        cfg_set_double(opts, "PhiMax", PhiMax);
        abn_write(fcells, &cells[0], Ncells, "2i8d", opts);
        cfg_destroy(opts);
    }
    else if(coordsys == "cartesian") {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz,XMin,XMax,YMin,YMax,ZMin,ZMax", ",")) {
            fprintf(stderr, "basis: must provide config options N{x,y,z} and {X,Y,Phi}{Z,Max}\n");
//            return 1;
            opsec_exit(1);
        }

        int Nx = cfg_get_int(cfg, "Nx");
        int Ny = cfg_get_int(cfg, "Ny");
        int Nz = cfg_get_int(cfg, "Nz");
        double XMin = cfg_get_double(cfg, "XMin");
        double XMax = cfg_get_double(cfg, "XMax");
        double YMin = cfg_get_double(cfg, "YMin");
        double YMax = cfg_get_double(cfg, "YMax");
        double ZMin = cfg_get_double(cfg, "ZMin");
        double ZMax = cfg_get_double(cfg, "ZMax");

        /* Populate list of cells */
        Cell c;
        int d, e, f;
        #pragma omp parallel for private(c,d,e,f)
        for(d = 0; d < Nx; d++) {
            c.xmin = XMin +     d*(XMax - XMin)/Nx;   // this choice of linear radial spacing is consistent with the assumptions in dot.cpp
            c.xmax = XMin + (d+1)*(XMax - XMin)/Nx;
            for(e = 0; e < Ny; e++) {
                c.ymin = YMin +     e*(YMax - YMin)/Ny;
                c.ymax = YMin + (e+1)*(YMax - YMin)/Ny;
                for(f = 0; f < Nz; f++) {
                    c.zmin = ZMin +     f*(ZMax - ZMin)/Nz;
                    c.zmax = ZMin + (f+1)*(ZMax - ZMin)/Nz;
                    c.G = (d*Ny + e)*Nz + f;

                    if(FinalizeCellC(c, nbar)) {
                        #pragma omp critical (cells_update)
                        {
                            c.a = Ncells++;
                            cells.push_back(c);
                        }
                    }
                }
            }
        }

        printf("Got %d non-empty cells out of %d possible.\n", Ncells, Nx*Ny*Nz);
        printf("Writing list of cells to '%s'.\n", cellfile);

        fprintf(fcells, "# struct Cell {\n");
        fprintf(fcells, "#     int a;\n");
        fprintf(fcells, "#     int G;\n");
        fprintf(fcells, "#     double xmin, xmax;\n");
        fprintf(fcells, "#     double ymin, ymax;\n");
        fprintf(fcells, "#     double zmin, zmax;\n");
        fprintf(fcells, "#     double Veff;\n");
        fprintf(fcells, "#     double Nbar;\n");
        fprintf(fcells, "# };\n");

        Config opts = cfg_new();
        cfg_set(opts, "coordsys", "cartesian");
        cfg_set_int(opts, "Ncells", Ncells);
        cfg_set_int(opts, "Nx", Nx);
        cfg_set_int(opts, "Ny", Ny);
        cfg_set_int(opts, "Nz", Nz);
        cfg_set_double(opts, "XMin", XMin);
        cfg_set_double(opts, "XMax", XMax);
        cfg_set_double(opts, "YMin", YMin);
        cfg_set_double(opts, "YMax", YMax);
        cfg_set_double(opts, "ZMin", ZMin);
        cfg_set_double(opts, "ZMax", ZMax);
        abn_write(fcells, &cells[0], Ncells, "2i8d", opts);
        cfg_destroy(opts);
    }
    
    return 0;
}
