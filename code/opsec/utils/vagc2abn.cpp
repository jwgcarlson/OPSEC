/* vagc2abn
 *
 * Takes a lss.$sample$letter$post.dat file of galaxies from the NYU VAGC and
 * converts it to a .abn galaxy file to be ingested by OPSEC.  The NYU VAGC
 * galaxy files are ASCII files organized into the following columns of data:
 *  (FINISH ME)
 *   . . . . . . . . . . .
 *
 * Uses WMAP7+BAO parameters for the redshift-to-distance conversion. */

#include <cmath>
#include <cstdio>
#include <unistd.h>

#include <vector>
using std::vector;

#include "abn.h"
#include "cfg.h"
#include "Survey.h"

const char* usage =
    "Usage: vagc2abn <.dat file> <.abn file>\n";

const double OmegaM = 0.272;
const double OmegaL = 0.728;
const double OmegaK = 1 - OmegaM - OmegaL;

double E(double z) {
    return sqrt(OmegaM*pow(1+z, 3) + OmegaK*pow(1+z, 2) + OmegaL);
}

/* Convert redshift cz (in km/s) to comoving distance \chi (in Mpc/h). */
double cz2chi(double cz) {
    double z = cz/299792.458;   // Upper integration limit
    double zp = 0.;             // Integration variable z'
    double chi = 0.;            // \chi(z)
    double dz0 = 0.001;         // Default integration step
    double dz;
    while(zp < z) {
        double dz = fmin(dz0, z - zp);
        chi += dz*(0.5/E(zp) + 0.5/E(zp+dz));   // Trapezoid rule
        zp += dz0;
    }
    chi *= 2997.92458;          // Convert to Mpc/h
    return chi;
}

/* Construct galaxy list */
void read_galaxies(FILE* fdat, vector<Galaxy>& galaxies) {
    int indx, sector, mregion;
    float ra, dec, cz, fgotten, selectfn;
    Galaxy g;
    while(fscanf(fdat, "%d %d %d %f %f %f %f %f", &indx, &sector, &mregion, &ra, &dec, &cz, &fgotten, &selectfn) == 8) {
        g.r = cz2chi(cz);
        g.mu = cos(M_PI/180. * (90. - dec));
        g.phi = M_PI/180. * ra;
        g.w = 1;        // ?
        galaxies.push_back(g);
    }
}

int main(int argc, char* argv[]) {
    /* Parse command line switches */
    int opt;
    const char* optstring = "";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        default:
            fputs(usage, stderr);
            return 1;
        }
    }

    if(argc != 3) {
        fputs(usage, stderr);
        return 1;
    }

    const char* datfile = argv[1];
    const char* abnfile = argv[2];

    vector<Galaxy> galaxies;
    FILE* fdat = fopen(datfile, "r");
    read_galaxies(fdat, galaxies);

    /* Write galaxies to file */
    FILE* fout = fopen(abnfile, "w");
    fprintf(fout, "# Galaxies converted to .abn format from '%s'\n", datfile);
    Config opts = cfg_new();
    cfg_set_int(opts, "ngals", (int) galaxies.size());
    abn_write(fout, &galaxies[0], galaxies.size(), "4f", opts);
    cfg_destroy(opts);
    fclose(fout);

    return 0;
}
