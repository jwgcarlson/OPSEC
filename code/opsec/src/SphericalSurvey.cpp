#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vector>
using std::vector;

#include "SphericalSurvey.h"
#include "Spline.h"
#include "abn.h"
#include "chealpix.h"


/* Selection function in spherical coordinates, using HEALPix map for angular
 * completeness map. */
struct SphericalSelectionFuncImpl : public SelectionFuncImpl {
    long nside;                 // HEALPix resolution parameter
    int nest;                   // pixel ordering: 1 for NESTED, 0 for RING
    float* mask;                // angular completeness at each pixel
    Spline radial_nbar;         // radial selection function

    SphericalSelectionFuncImpl(long nside_, int nest_, float* mask_, const Spline& radial_nbar_) {
        nside = nside_;
        nest = nest_;
        mask = mask_;
//        mask = (float*)malloc(12*nside*nside*sizeof(float));
//        memcpy(mask, mask_, 12*nside*nside*sizeof(float));
        radial_nbar = radial_nbar_;
    }

    ~SphericalSelectionFuncImpl() {
        free(mask);
    }

    double evaluate(double r, double mu, double phi) {
        double rad = radial_nbar(r);
        double theta = acos(mu);
        long ip;
        if(nest)
            ang2pix_nest(nside, theta, phi, &ip);
        else
            ang2pix_ring(nside, theta, phi, &ip);
        double ang = mask[ip];
        return ang*rad;
    }
};


SphericalSurvey::SphericalSurvey(Config cfg) {
    galfile = cfg_get(cfg, "galfile");
    maskfile = cfg_get(cfg, "maskfile");
    radial = cfg_get(cfg, "radial");
}

SphericalSurvey::~SphericalSurvey() {
}

SelectionFunc SphericalSurvey::GetSelectionFunction() {
    /* Read HEALPix file */
    long nside;
    char coordsys[9];
    char ordering[32];
    int nest;
    float* mask = read_healpix_map(maskfile.c_str(), &nside, coordsys, ordering);
    if(mask == NULL) {
        fprintf(stderr, "SphericalSurvey: could not read HEALPix file '%s'\n", maskfile.c_str());
        return NULL;
    }
    else {
        if(strncmp(ordering, "NESTED", 6) == 0)
            nest = 1;
        else if(strncmp(ordering, "RING", 4) == 0)
            nest = 0;
        else {
            fprintf(stderr, "SphericalSurvey: unrecognized HEALPix ordering: %s\n", ordering);
            return NULL;
        }
    }

    Spline radial_nbar;
    if(radial == "compute") {
#if 0
        double A = ComputeWeightedArea(nside, mask);
        printf("Weighted sky area coverage = %g steradians\n", A);

        /* Compute radial selection function from galaxies. */
        int nbins = 20;
        vector<double> rbins(nbins);
#endif
        fprintf(stderr, "SphericalSurvey: 'radial = compute' not implemented yet (FIXME)\n");
        return NULL;
    }
    else {
        /* Read radial selection function from file */
        radial_nbar = LinearSpline(radial.c_str());
        return new SphericalSelectionFuncImpl(nside, nest, mask, radial_nbar);
    }
}

Galaxy* SphericalSurvey::GetGalaxies(int* ngals) {
    FILE* fgals = fopen(galfile.c_str(), "r");
    if(fgals == NULL) {
        fprintf(stderr, "SphericalSurvey: could not open file '%s'\n", galfile.c_str());
        return NULL;
    }

    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    if(abn_read_header(fgals, &n, &size, &endian, fmt, NULL) != 0
       || size != sizeof(Galaxy)
       || strcmp(fmt, "4f") != 0)
    {
        fprintf(stderr, "SphericalSurvey: error reading galaxies from '%s'\n", galfile.c_str());
        fclose(fgals);
        return NULL;
    }

    Galaxy* gals = (Galaxy*) malloc(n*sizeof(Galaxy));
    if(!gals) {
        fprintf(stderr, "BoxSurvey: could not allocate memory for %zd galaxies\n", n);
        fclose(fgals);
        return NULL;
    }

    size_t nread = fread((void*) gals, size, n, fgals);
    if(nread != n) {
        fprintf(stderr, "SphericalSurvey: expecting %zd galaxies from '%s', got %zd\n", n, galfile.c_str(), nread);
        fclose(fgals);
        return NULL;
    }

    fclose(fgals);
    if(ngals)
        *ngals = (int)n;
    return gals;
}

double SphericalSurvey::ComputeWeightedArea(long nside, float* mask) {
    assert(mask != NULL);

    long npix = 12*nside*nside;
    double A = 0.;
    for(long i = 0; i < npix; i++)
        A += mask[i];
    A *= 4*M_PI/npix;
    return A;
}

void SphericalSurvey::ComputeRadialNbar(
    int ngals,
    const Galaxy* gals,
    double A,               // completeness-weighted sky area in steradians
    int nr,
    const double* r,
    double* nbar)
{
    /* (Weighted) number of galaxies in each radial bin */
    double* counts = (double*) calloc(nr-1, sizeof(double));

    for(int g = 0; g < ngals; g++) {
        if(gals[g].r < r[0])
            continue;
        for(int i = 0; i < nr-1; i++) {
            if(gals[g].r < r[i+1]) {
                counts[i] += gals[g].w;
                break;
            }
        }
    }

    for(int i = 0; i < nr-1; i++) {
        double Vshell = (pow(r[i+1], 3) - pow(r[i], 3))/3. * A;
        nbar[i] = counts[i] / Vshell;
    }
    free(counts);
}
