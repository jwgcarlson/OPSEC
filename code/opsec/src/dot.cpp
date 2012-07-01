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
#include "Model.h"
#include "SeparationFunc.h"
#include "Spline.h"
#include "Survey.h"
#include "XiFunc.h"
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
    "  pixfile=FILE  Write pixel values to FILE (required)\n"
    "  cellfile=FILE Read list of cells from FILE (required)\n"
    "  modefile=FILE Read KL mode coefficients from FILE (required)\n"
    "  Nmodes=NUM    Number of KL modes\n";

int main(int argc, char* argv[]) {
    Config cfg = opsec_init(argc, argv, usage);

    /* Read configuration options */
    if(cfg_missing_keys(cfg, "countfile,cellfile,modefile,pixfile,Nmodes")) {
        opsec_error("dot: missing configuration options\n");
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
        opsec_error("dot: could not open '%s' for writing\n", pixfile);
        opsec_exit(1);
    }

    /* Read basis cells */
    int Ncells;
    Cell* cells = ReadCells(cellfile, &Ncells);
    if(cells == NULL) {
        opsec_error("dot: error reading cells\n");
        opsec_exit(1);
    }

    /* Read KL modes (expressed in terms of basis cells) */
    int Nmodes_check, Ncells_check;
    FILE* fmodes = ReadModesHeader(modefile, &Nmodes_check, &Ncells_check);
    if(Nmodes != Nmodes_check || Ncells != Ncells_check) {
        fprintf(stderr, "dot: inconsistent modes file\n");
        opsec_exit(1);
    }

    opsec_debug("Read %d modes from '%s'\n", Nmodes, cellfile);

    /* Load survey */
    Survey* survey = InitializeSurvey(cfg);
    if(!survey)
        opsec_exit(1);

    /* Read in galaxies */
    vector<Galaxy> gals;
    survey->GetGalaxies(gals);
    int Ngals = (int) gals.size();
    if(Ngals == 0) {
        opsec_error("dot: error reading galaxies\n");
        opsec_exit(1);
    }
    opsec_debug("Read %d galaxies\n", Ngals);

    /* Construct a mapping between the grid index e = (d1*N2 + d2)*N3 + d3 of a
     * cell and its index (a) within the list of non-empty cells. */
    int N1 = survey->N1;
    int N2 = survey->N2;
    int N3 = survey->N3;
    vector<int> cellmap(N1*N2*N3, -1);
    for(int a = 0; a < Ncells; a++) {
        const Cell& c = cells[a];
        opsec_assert(a == c.a);
        int e = (c.d1*N2 + c.d2)*N3 + c.d3;
        cellmap[e] = a;
    }

    /* Number of galaxies found in each cell. The last element (counts[Ncells])
     * gives the number of galaxies that don't lie in any cell.
     * Note that each galaxy carries a real-valued weight, so the count for any
     * given cell need not be integral. */
    double* counts = (double*) opsec_calloc(Ncells+1, sizeof(double));

    /* Count the (weighted) number of galaxies within each non-empty cell */
    for(vector<Galaxy>::const_iterator g = gals.begin(); g != gals.end(); g++) {
        int e = survey->GetGridIndex(g->x1, g->x2, g->x3);
        int a = (e > 0) ? cellmap[e] : Ncells;
        counts[a] += g->w;
    }

    /* x_a = \int d^3x \phi_a(\vec{x}) [n(\vec{x}) - \bar{n}(\vec{x})]
     *     = \bar{N}_a^{-1/2} [N_a - \bar{N}_a] */
    vector<double> x(Ncells, 0);
    for(int a = 0; a < Ncells; a++)
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

#if 0
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
    double S = (coordsys == "spherical") ? ComputeSignalS(cells[0], cells[0], xi, survey->GetSeparationFunction())
                                         : ComputeSignalC(cells[0], cells[0], xi, survey->GetSeparationFunction());
    for(int a = 0; a < Ncells; a++) {
//        double S = (coordsys == "spherical") ? ComputeSignalS(cells[a], cells[a], xi)
//                                             : ComputeSignalC(cells[a], cells[a], xi);
        fprintf(fcounts, "%6d %8.4f %8.4f %8.4f\n", a, counts[a], cells[a].Nbar, cells[a].Nbar*(1 + cells[a].Nbar*S));
    }
#endif

    /* Clean up */
    free(cells);
    free(counts);
    return 0;
}
