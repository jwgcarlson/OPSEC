#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vector>
using std::vector;

#include "SphericalSurvey.h"
#include "SelectionFunc.h"
#include "SeparationFunc.h"
#include "Spline.h"
#include "abn.h"
#include "chealpix.h"


/* Separation between two points in spherical coordinates. */
struct SphericalSeparationFuncImpl : public SeparationFuncImpl {
    double r(const Point& p1, const Point& p2) {
        double cosgamma = sqrt((1-p1.mu*p1.mu) * (1-p2.mu*p2.mu)) * cos(p1.phi - p2.phi) + p1.mu*p2.mu;
        return sqrt(fmax(p1.r*p1.r + p2.r*p2.r - 2*p1.r*p2.r*cosgamma, 0.0));
    }

    void rmu(const Point& p1, const Point& p2, double& r, double& mu) {
        double cosgamma = sqrt((1-p1.mu*p1.mu) * (1-p2.mu*p2.mu)) * cos(p1.phi - p2.phi) + p1.mu*p2.mu;
        r = sqrt(fmax(p1.r*p1.r + p2.r*p2.r - 2*p1.r*p2.r*cosgamma, 0.0));
        if(r == 0)
            mu = 0;
        else {
            double singamma = sqrt(1 - fmax(1, cosgamma*cosgamma));
            double sinbeta = singamma * p2.r/r;         // \sin\gamma / r = \sin\beta / b
            double cosbeta = sqrt(1 - fmax(1, sinbeta*sinbeta));
            double cosgamma2 = sqrt((1 + cosgamma)/2.); // \cos(\gamma/2) = \sqrt{(1 + \cos\gamma)/2}
            double singamma2 = sqrt((1 - cosgamma)/2.); // \sin(\gamma/2) = \sqrt{(1 - \cos\gamma)/2}
            mu = cosgamma2*cosbeta - singamma2*sinbeta; // \mu = \cos\phi = \cos(\gamma/2 + \beta)
        }
    }

    void rab(const Point& p1, const Point& p2, double& r, double& a, double& b) {
        double cosgamma = sqrt((1-p1.mu*p1.mu) * (1-p2.mu*p2.mu)) * cos(p1.phi - p2.phi) + p1.mu*p2.mu;
        r = sqrt(fmax(p1.r*p1.r + p2.r*p2.r - 2*p1.r*p2.r*cosgamma, 0.0));
        a = p1.r;
        b = p2.r;
    }
};


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


SphericalSurvey::SphericalSurvey(Config cfg)
    : Survey(CoordSysSpherical, cfg)
{
    galfile = cfg_get(cfg, "galfile");
    maskfile = cfg_get(cfg, "maskfile");
    radial = cfg_get(cfg, "radial");
}

SphericalSurvey::~SphericalSurvey() {
}


SeparationFunc SphericalSurvey::GetSeparationFunction() {
    return SeparationFunc(new SphericalSeparationFuncImpl());
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
        return SelectionFunc(new SphericalSelectionFuncImpl(nside, nest, mask, radial_nbar));
    }
}

void SphericalSurvey::GetGalaxies(std::vector<Galaxy>& gals) {
    size_t norig = gals.size();

    FILE* fgals = fopen(galfile.c_str(), "r");
    if(fgals == NULL) {
        fprintf(stderr, "SphericalSurvey: could not open file '%s'\n", galfile.c_str());
        return;
    }

    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    if(abn_read_header(fgals, &n, &size, &endian, fmt, NULL) != 0
       || size != sizeof(Galaxy)
       || strcmp(fmt, "4f") != 0)
    {
        fprintf(stderr, "SphericalSurvey: error reading galaxies from '%s'\n", galfile.c_str());
        fclose(fgals);
        return;
    }

    /* Make room for n additional Galaxy objects */
    gals.resize(norig + n);

    size_t nread = fread((void*) &gals[norig], size, n, fgals);
    if(nread != n)
        fprintf(stderr, "SphericalSurvey: expecting %zd galaxies from '%s', got %zd\n", n, galfile.c_str(), nread);

    fclose(fgals);
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
