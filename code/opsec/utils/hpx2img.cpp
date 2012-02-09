/* hpx2img
 *
 * Converts a HEALPix completeness map into an image file.
 *
 * TODO: support multiple projections
 */

#include <cmath>
#include <cstdio>
#include <unistd.h>

#include "chealpix.h"

#define cimg_display 0
#define cimg_use_png
#include "CImg.h"
using namespace cimg_library;


/* Compute the Hammer-Aitoff projection (x,y) of the spherical coordinates
 * (mu,phi).  On return, x and y lie in [-1,+1]. */
void hammer(double mu, double phi, double* x, double* y) {
    mu = fmin(+1, fmax(-1, mu));        // clamp mu to [-1,+1]
    double sintheta = sqrt(1 - mu*mu);  // sin(theta) = cos(latitude)
    double longitude = phi - M_PI;
    double z = sqrt(1 + sintheta*cos(longitude/2));
    *x = sintheta*sin(longitude/2)/z;
    *y = mu/z;
}

/* Compute the inverse Hammer-Aitoff projection (mu,phi) of the map coordinates
 * (x,y), normalized to lie in [-1,+1].  If the map coordinates lie outside the
 * projection, mu and phi are untouched and a non-zero value is returned.
 * Otherwise 0 is returned, mu lies in [-1,+1] and phi lies in [0,2*pi). */
int inverse_hammer(double x, double y, double* mu, double* phi) {
    if(1 - x*x - y*y <= 0)
        return 1;
    double z = sqrt(1 - x*x/2. - y*y/2.);
    double longitude = 2*atan(M_SQRT2*x*z/(2*z*z-1));
    double latitude = asin(M_SQRT2*y*z);
    *mu = sin(latitude);
    *phi = longitude + M_PI;
    return 0;
}

int main(int argc, char* argv[]) {
    const char* hpxfile = NULL;
    const char* imgfile = NULL;
    int width = 800;
    int height = 400;
    int xmargin = 0;
    int ymargin = 0;
    int nsamp = 4;      // number of times to sample each pixel (in each direction)
    int raw = 0;        // if nonzero, assume hpxfile contains raw HEALPix data with nside=raw

    /* Parse command line switches */
    int opt;
    while((opt = getopt(argc, argv, "hg:m:s:r:")) != -1) {
        switch(opt) {
        case 'h':
            printf("Usage: %s [-g width,height] [-m xmargin,ymargin] [-s nsamp] hpxfile imgfile\n", argv[0]);
            exit(0);
            break;
        case 'g':
            if(sscanf(optarg, "%d,%d", &width, &height) != 2) {
                fprintf(stderr, "Invalid argument to -g\n");
                exit(1);
            }
            break;
        case 'm':
            if(sscanf(optarg, "%d,%d", &xmargin, &ymargin) != 2) {
                fprintf(stderr, "Invalid argument to -m\n");
                exit(1);
            }
            break;
        case 's':
            if(sscanf(optarg, "%d", &nsamp) != 1) {
                fprintf(stderr, "Invalid argument to -s\n");
                exit(1);
            }
            break;
        case 'r':
            if(sscanf(optarg, "%d", &raw) != 1) {
                fprintf(stderr, "Invalid argument to -r\n");
                exit(1);
            }
            break;
        default:
            fprintf(stderr, "Illegal argument -%c, try `%s -h`\n", opt, argv[0]);
            exit(1);
            break;
        }
    }

    /* Parse additional command line options */
    if(argc - optind != 2) {
        fprintf(stderr, "Usage error, try `%s -h`\n", argv[0]);
        exit(1);
    }
    hpxfile = argv[optind];
    imgfile = argv[optind+1];

    long nside, npix, ipix;
    float* map;

    /* If the string hpxfile is just a number in [0,1], assume a constant completeness map */
    float v;
    if(sscanf(hpxfile, "%f", &v) == 1 && v >= 0. && v <= 1.) {
        nside = 1;
        npix = 12;
        map = (float*) malloc(npix*sizeof(float));
        for(ipix = 0; ipix < npix; ipix++)
            map[ipix] = v;
    }
    else {
        printf("Reading HEALPix map from %s\n", hpxfile);
        if(raw == 0) {
            /* Read from HEALPix data from standard FITS file */
            char coordsys[9], ordering[9];
            map = read_healpix_map(hpxfile, &nside, coordsys, ordering);
            if(map == NULL) {
                fprintf(stderr, "could not read HEALPix map from %s\n", hpxfile);
                return 1;
            }
            npix = 12*nside*nside;
        }
        else {
            nside = raw;
            npix = 12*nside*nside;
            map = (float*) malloc(npix*sizeof(float));
            FILE* fhpx = fopen(hpxfile, "r");
            if(fhpx == NULL) {
                fprintf(stderr, "could not open %s\n", hpxfile);
                return 1;
            }
            if(fread(map, sizeof(float), npix, fhpx) != (size_t)npix) {
                fprintf(stderr, "error reading raw HEALPix data from %s\n", hpxfile);
                return 1;
            }
            fclose(fhpx);
        }
    }

    /* Save raw HEALPix data as a proper HEALPix/FITS file */
//    float* fmap = (float*) malloc(npix*sizeof(float));
//    for(long i = 0; i < npix; i++)
//        fmap[i] = map[i];
//    write_healpix_map(fmap, nside, "goodmask.hpx", 1, "G");

#if 0
    /* Normalize map */
    float pixmin = 1e10, pixmax = -1e10;
    for(ipix = 0; ipix < npix; ipix++) {
        if(map[ipix] < pixmin) pixmin = map[ipix];
        if(map[ipix] > pixmax) pixmax = map[ipix];
    }
#endif

    /* Initialize image to transparent */
    CImg<unsigned char> image(width, height, 1, 4, 0);

    printf("Painting completeness map\n");

    /* Color each pixel of the image according to map values */
    double x, y, mu, theta, phi;
    for(int i = xmargin; i < width - xmargin; i++) {
        for(int j = ymargin; j < height - ymargin; j++) {
            /* Sample each pixel multiple times */
            double value = 0; // average completeness value within pixel
            int ninside = 0;  // number of sample points within the Hammer-Aitoff projection region
            for(int is = 0; is < nsamp; is++) {
                for(int js = 0; js < nsamp; js++) {
                    x = -1. + 2.*(i + (is+0.5)/nsamp - xmargin)/(width - 2.*xmargin);
                    y =  1. - 2.*(j + (js+0.5)/nsamp - ymargin)/(height - 2.*ymargin);  // flip, since pixels run top-to-bottom
                    if(inverse_hammer(x, y, &mu, &phi) == 0) {
                        theta = acos(mu);
                        ang2pix_nest(nside, theta, phi, &ipix);
                        value += map[ipix];
                        ninside++;
                    }
                }
            }
            value /= ninside;

            if(ninside > 0) {
                value = fmax(0, fmin(1, value)); // clamp to [0,1]
                /* Convert pixel value to a color */
                image(i,j,0,0) = 255*(1 - value);
                image(i,j,0,1) = 255*(1 - value);
                image(i,j,0,2) = 255*(1 - value);
                image(i,j,0,3) = 255;   // opaque
            }
        }
    }

#if 0
    /* Draw elliptical outline around map */
    int npoints = 1024;
    CImg<double> line(npoints, 2);
    line(0,0) = width/2.;
    line(0,1) = height-ymargin;
    for(int i = 0; i < npoints/2; i++) {
        mu = 1. - i*4./npoints;
        hammer(mu, 2*M_PI, &x, &y);
        line(i,0) = xmargin + (width-2*xmargin)*(x + 1.)/2.;
        line(i,1) = ymargin + (height-2*ymargin)*(y + 1.)/2.;
    }
    line(npoints/2,0) = width/2.;
    line(npoints/2,0) = ymargin;
    for(int i = 0; i < npoints/2; i++) {
        mu = -1. + i*4./npoints;
        hammer(mu, 0, &x, &y);
        line(i+npoints/2,0) = xmargin + (width-2*xmargin)*(x + 1.)/2.;
        line(i+npoints/2,1) = ymargin + (height-2*ymargin)*(y + 1.)/2.;
    }
    image.draw_polygon(line, "black", 1.0, ~0U);
#endif

    /* Write image to file */
    image.save(imgfile);

    free(map);
    return 0;
}
