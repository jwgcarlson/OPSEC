/* makeradial
 *
 * Compute the radial number density $\bar{n}(r)$ for a given survey. */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vector>
using std::vector;

#include "SphericalSurvey.h"
#include "abn.h"
#include "cfg.h"
#include "chealpix.h"

int main(int argc, char* argv[]) {
    if(argc != 4) {
        fprintf(stderr, "Usage: makeradial <galfile> <maskfile> <radialfile>\n");
        return 1;
    }

    const char* galfile = argv[1];
    const char* maskfile = argv[2];
    const char* radialfile = argv[3];

    /* Load mask */
    long nside;
    char coordsys[9], ordering[32];
    float* mask = read_healpix_map(maskfile, &nside, coordsys, ordering);
    int nest = strncmp(ordering, "NESTED", 6) == 0;

    /* Load galaxies */
    FILE* fgals = fopen(galfile, "r");
    if(fgals == NULL) {
        fprintf(stderr, "makeradial: could not open %s\n", galfile);
        return 1;
    }
    Galaxy* gals = NULL;
    size_t ngals;
    if(abn_read(fgals, (void**) &gals, &ngals, NULL, NULL, NULL, NULL) != 0) {
        fprintf(stderr, "makeradial: error reading galaxies from %s\n", galfile);
        return 1;
    }
    fclose(fgals);

    /* Compute sky area */
    double A = SphericalSurvey::ComputeWeightedArea(nside, mask);
    printf("Effective sky area = %g steradians\n", A);

    /* Find maximum r, and round to a multiple of 100 Mpc/h */
    double rmax = 0.;
    for(int g = 0; g < ngals; g++)
        rmax = (gals[g].r > rmax) ? gals[g].r : rmax;
    rmax = 100*(int(rmax/100) + 1);
    int nr = 1 + 2*(int(rmax/100) + 1); // 0, 50, 100, 150, ..., rmax-50, rmax

    /* Compute radial number density */
    vector<double> rr(nr), nbar(nr);
    for(int i = 0; i < nr; i++)
        rr[i] = 50.*i;
    SphericalSurvey::ComputeRadialNbar((int) ngals, gals, A, nr, &rr[0], &nbar[0]);

    /* Write to file */
    FILE* fradial = fopen(radialfile, "w");
    for(int i = 0; i < nr; i++)
        fprintf(fradial, "%g %g\n", rr[i], nbar[i]);
    fclose(fradial);

    printf("Wrote radial number density to %s\n", radialfile);

    free(mask);
    return 0;
}
