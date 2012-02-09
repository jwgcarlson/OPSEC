/* dot
 * 
 * Compute pixel values $y_i$ by dotting the observed galaxy density field with
 * the KL mode functions.
 *
 * This program is designed to be run on a single node, using OpenMP for
 * shared-memory parallelism. */

/* TODO:
 * - make cellmap use longs (or a guaranteed 64-bit integer), in case cell grid
 *   is really large */

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

#include "Cell.h"
#include "Spline.h"
#include "Survey.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"
#include "sig.h"

const char* usage =
    "Usage: dot [SWITCHES] [OPTIONS]\n"
    "Switches:\n"
    "  -h            Display this help\n"
    "  -c FILE       Read additional configuration options from FILE\n"
    "Configuration options:\n"
    "  coordsys=TYPE Coordinate system; either spherical or cartesian (required)\n"
    "  survey=NAME   Survey configuration (required)\n"
    "  pixfile=FILE  Write pixel values to FILE (required)\n"
    "  cellfile=FILE Read list of cells from FILE (required)\n"
    "  modefile=FILE Read KL mode coefficients from FILE (required)\n"
    "  Nmodes=NUM    Number of KL modes\n";

int main(int argc, char* argv[]) {
    Config cfg = opsec_init(argc, argv, usage);

    /* Read configuration options */
    if(!cfg_has_keys(cfg, "countfile,cellfile,modefile,pixfile,Nmodes", ",")) {
        fprintf(stderr, "dot: missing configuration options\n");
        fputs(usage, stderr);
        opsec_exit(1);
    }
    const char* countfile = cfg_get(cfg, "countfile");
    const char* cellfile = cfg_get(cfg, "cellfile");
    const char* modefile = cfg_get(cfg, "modefile");
    const char* pixfile = cfg_get(cfg, "pixfile");
    int Nmodes = cfg_get_int(cfg, "Nmodes");

    /* Make sure we can write to the output file */
    FILE* fpix = fopen(pixfile, "w");
    if(fpix == NULL) {
        fprintf(stderr, "dot: could not open '%s' for writing\n", pixfile);
        opsec_exit(1);
    }

    /* Read basis cells */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL) {
        fprintf(stderr, "dot: error reading cells\n");
        opsec_exit(1);
    }

    /* Read KL modes (expressed in terms of basis cells) */
    int Nmodes_check, Ncells_check;
    FILE* fmodes = ReadModesHeader(modefile, &Nmodes_check, &Ncells_check);
    if(Nmodes != Nmodes_check || Ncells != Ncells_check) {
        fprintf(stderr, "dot: inconsistent modes file\n");
        opsec_exit(1);
    }

    printf("(TRACE) Read %d modes from '%s'\n", Nmodes, cellfile); fflush(stdout);

    /* Read galaxies */
    Survey* survey = InitializeSurvey(cfg);
    if(survey == NULL) {
        fprintf(stderr, "dot: error initializing survey\n");
        opsec_exit(1);
    }

    printf("Initialized survey\n");
    fflush(stdout);

    int Ngals = 0;
    Galaxy* gals = survey->GetGalaxies(&Ngals);
    printf("Ngals = %d\n", Ngals);
    if(gals == NULL || Ngals == 0) {
        fprintf(stderr, "dot: error reading galaxies\n");
        opsec_exit(1);
    }

    printf("Read %d galaxies\n", Ngals);
    fflush(stdout);

    /* Number of galaxies found in each cell. The last element (counts[Ncells])
     * gives the number of galaxies that don't lie in any cell.
     * Note that each galaxy carries a real-valued weight, so the count for any
     * given cell need not be integral. */
    double* counts = (double*)opsec_calloc(Ncells+1, sizeof(double));

    string coordsys(cfg_get(cfg, "coordsys"));
    if(coordsys != "spherical" && coordsys != "cartesian") {
        fprintf(stderr, "dot: missing or invalid config option: coordsys = %s\n", coordsys.c_str());
        opsec_exit(1);
    }

    int a;
    if(coordsys == "spherical") {
        if(!cfg_has_keys(cfg, "Nr,Nmu,Nphi,RMin,RMax,MuMin,MuMax,PhiMin,PhiMax", ",")) {
            fprintf(stderr, "dot: must provide config options N{r,mu,phi} and {R,Mu,Phi}{Min,Max}\n");
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

        /* Create mapping between cell grid location G and cell number a, with
         * 0 <= a < Ncells corresponding to actual cells.  A value of a = Ncells
         * means the cell was ignored (i.e. selection function vanishes over entire
         * cell volume). */
        vector<int> cellmap(Nr*Nmu*Nphi, Ncells);
        for(a = 0; a < Ncells; a++) {
            assert(a == cells[a].a);
            cellmap[cells[a].G] = a;
        }

        /* Count number of galaxies in each cell */
        double r, mu, phi, w;
        int d, e, f, G, g;
//        #pragma omp parallel for private(g,r,mu,phi,d,e,f,G,a)
        for(g = 0; g < Ngals; g++) {
            r = gals[g].r;
            mu = gals[g].mu;
            phi = gals[g].phi;
            w = gals[g].w;

            if((RMin <= r && r < RMax) && (MuMin <= mu && mu < MuMax) && (PhiMin <= phi && phi < PhiMax)) {
                /* Galaxy lies within cell grid.  Find which cell it lies in. */
                d = (int) floor(Nr * (r - RMin)/(RMax - RMin));
                e = (int) floor(Nmu * (mu - MuMin)/(MuMax - MuMin));
                f = (int) floor(Nphi * (phi - PhiMin)/(PhiMax - PhiMin));
                G = (d*Nmu + e)*Nphi + f;
                a = cellmap[G];
//                #pragma omp atomic
                counts[a] += w;
            }
            else {
                /* Cell lies outside grid.  Count this as lying in an empty cell. */
//                #pragma omp atomic
                counts[Ncells] += w;
            }
        }
    }
    else if(coordsys == "cartesian") {
        if(!cfg_has_keys(cfg, "Nx,Ny,Nz,XMin,XMax,YMin,YMax,ZMin,ZMax", ",")) {
            fprintf(stderr, "dot: must provide config options N{x,y,z} and {X,Y,Phi}{Z,Max}\n");
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

        /* Create mapping between cell grid location G and cell number a, with
         * 0 <= a < Ncells corresponding to actual cells.  A value of a = Ncells
         * means the cell was ignored (i.e. selection function vanishes over entire
         * cell volume). */
        vector<int> cellmap(Nx*Ny*Nz, Ncells);
        for(a = 0; a < Ncells; a++) {
            assert(a == cells[a].a);
            cellmap[cells[a].G] = a;
        }

        /* Count number of galaxies in each cell */
        double x, y, z, w;
        int d, e, f, G, g;
//        #pragma omp parallel for private(g,x,y,z,d,e,f,G,a)
        for(g = 0; g < Ngals; g++) {
            x = gals[g].x;
            y = gals[g].y;
            z = gals[g].z;
            w = gals[g].w;

            if((XMin <= x && x < XMax) && (YMin <= y && y < YMax) && (ZMin <= z && z < ZMax)) {
                /* Galaxy lies within cell grid.  Find which cell it lies in. */
                d = (int) floor(Nx * (x - XMin)/(XMax - XMin));
                e = (int) floor(Ny * (y - YMin)/(YMax - YMin));
                f = (int) floor(Nz * (z - ZMin)/(ZMax - ZMin));
                G = (d*Ny + e)*Nz + f;
                a = cellmap[G];
//                #pragma omp atomic
                counts[a] += w;
            }
            else {
                /* Cell lies outside grid.  Count this as lying in an empty cell. */
//                #pragma omp atomic
                counts[Ncells] += w;
            }
        }
    }

    /* x_a = \int d^3x \phi_a(\vec{x}) [n(\vec{x}) - \bar{n}(\vec{x})]
     *     = \bar{N}_a^{-1/2} [N_a - \bar{N}_a] */
    vector<double> x(Ncells, 0);
//    #pragma omp parallel for
    for(a = 0; a < Ncells; a++)
        x[a] = (counts[a] - cells[a].Nbar)/sqrt(cells[a].Nbar);

    /* y_i = \int d^3x \psi_i(\vec{x}) [n(\vec{x}) - \bar{n}(\vec{x})]
     *     = \sum_a B_{ia} \int d^3x \phi_a(\vec{x}) [n(\vec{x}) - \bar{n}(\vec{x})]
     *     = \sum_a B_{ia} \bar{N}_a^{-1/2} [N_a - \bar{N}_a] */
    vector<double> y(Nmodes, 0);

    /* Subtract mean and normalize to get pixel values. */
    vector<real> mode(Ncells, 0);
    int incx = 1, incy = 1;
    for(int i = 0; i < Nmodes; i++) {
        /* Read one mode at a time */
        size_t nread = fread(&mode[0], sizeof(real), Ncells, fmodes);
        assert(nread == Ncells);

        /* Compute y_i = B_{ia} x_a */
        y[i] = blas_dot(&Ncells, &mode[0], &incx, &x[0], &incy);
    }
    fclose(fmodes);

    /* Write pixel values to file */
    fprintf(fpix, "# Pixel values\n");
    abn_write(fpix, &y[0], Nmodes, "d", NULL);
    fclose(fpix);

    /* Load model */
    Model* model = InitializeModel(cfg);
    if(!model) {
        fprintf(stderr, "dot: could not load model\n");
        opsec_exit(1);
    }
    XiFunc xi = model->GetXi();

    /* Write counts to ASCII file */
    FILE* fcounts = fopen(countfile, "w");
    if(!fcounts) {
        fprintf(stderr, "dot: could not open '%s' for writing\n", countfile);
        opsec_exit(1);
    }
    fprintf(fcounts, "# Cell counts\n");
    fprintf(fcounts, "# %g (weighted) galaxies outside any cell\n", counts[Ncells]);
    fprintf(fcounts, "# Columns are:  cell number -- observed count -- expected count -- variance of count\n");
    double S = (coordsys == "spherical") ? ComputeSignalS(cells[0], cells[0], xi)
                                         : ComputeSignalC(cells[0], cells[0], xi);
    for(int a = 0; a < Ncells; a++) {
//        double S = (coordsys == "spherical") ? ComputeSignalS(cells[a], cells[a], xi)
//                                             : ComputeSignalC(cells[a], cells[a], xi);
        fprintf(fcounts, "%6d %8.4f %8.4f %8.4f\n", a, counts[a], cells[a].Nbar, cells[a].Nbar*(1 + cells[a].Nbar*S));
    }

    /* Clean up */
    free(cells);
    free(counts);
    free(gals);
    return 0;
}
