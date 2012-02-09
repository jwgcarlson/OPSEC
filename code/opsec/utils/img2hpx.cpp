/* img2hpx
 *
 * Converts an image file into a HEALPix completeness map.
 *
 * TODO: support multiple projections
 */

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <unistd.h>

#include "chealpix.h"

#define cimg_display 0
#define cimg_use_png
#include "CImg.h"
using namespace cimg_library;

enum ProjectionType {
    MuPhi,      // mu = 1 - 2*(y/h), phi = 2*pi*(x/w)
    ThetaPhi,   // theta = pi*(y/h), phi = 2*pi*(x/w)
    Aitoff,     // mu,phi = aitoff(2*(x/w)-1, 1-2*(y/h))
    Hammer      // mu,phi = hammer(2*(x/w)-1, 1-2*(y/h))
};

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

template<typename T>
bool ispow2(T x) {
    return !(x & (x-1));
}

void exit_msg(int status, FILE* stream, const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stream, format, ap);
    va_end(ap);
    exit(status);
}

int main(int argc, char* argv[]) {
    const char* imgfile = NULL;
    const char* hpxfile = NULL;
    int xmargin = 0;
    int ymargin = 0;
    long nside = 1024;
    ProjectionType projection;
    int raw = 0;        // if nonzero, assume write hpxfile as raw HEALPix data

    /* Parse command line switches */
    int opt;
    while((opt = getopt(argc, argv, "hm:n:p:r")) != -1) {
        switch(opt) {
        case 'h':
            exit_msg(0, stdout, "Usage: %s [-n nside] [-m xmargin,ymargin] [-p projection] imgfile hpxfile\n", argv[0]);
            break;
        case 'm':
            if(sscanf(optarg, "%d,%d", &xmargin, &ymargin) != 2)
                exit_msg(1, stderr, "Invalid argument to -m: %s\n", optarg);
            break;
        case 'n':
            if(sscanf(optarg, "%ld", &nside) != 1)
                exit_msg(1, stderr, "Invalid argument to -n: %s\n", optarg);
            break;
        case 'p':
            if(strcmp(optarg, "MuPhi") == 0)
                projection = MuPhi;
            else if(strcmp(optarg, "ThetaPhi") == 0)
                projection = ThetaPhi;
//            else if(strcmp(optarg, "Aitoff") == 0)
//                projection = Aitoff;
            else if(strcmp(optarg, "Hammer") == 0)
                projection = Hammer;
            else
                exit_msg(1, stderr, "Invalid projection: %s\n", optarg);
            break;
        case 'r':
            raw = 1;
            break;
        default:
            exit_msg(1, stderr, "Illegal argument -%c, try `%s -h`\n", opt, argv[0]);
            break;
        }
    }

    /* Parse additional command line options */
    if(argc - optind != 2)
        exit_msg(1, stderr, "Usage error, try `%s -h`\n", argv[0]);
    imgfile = argv[optind];
    hpxfile = argv[optind+1];

    /* Read in image file */
    CImg<unsigned char> image(imgfile);
    int width = image.width();
    int height = image.height();

    int w = width - 2*xmargin;
    int h = height - 2*ymargin;

    /* Sanity check arguments */
    if(nside < 1 || nside > 8192 || !ispow2(nside))
        exit_msg(1, stderr, "Error: nside must be a power of 2 between 1 and 8192 (not %ld)\n", nside);
    if(2*xmargin >= w || 2*ymargin >= h)
        exit_msg(1, stderr, "Error: margins are %dx%d for a %dx%d image\n", xmargin, ymargin, w, h);

    /* Initialize completeness mask */
    long ipix, npix = 12*nside*nside;
    float* map = (float*) calloc(npix, sizeof(float));
    int* count = (int*) calloc(npix, sizeof(int));

    /* Determine a sufficient sampling rate that each HEALPix pixel will be sampled at least once */
//    int nsamp = (int) ceil((4.*npix)/(1.*w*h));
//    printf("nsamp = %d\n", nsamp);
    int nsamp = 4;

    printf("Computing HEALPix completeness map\n");

    /* For each sampling of each pixel, update the corresponding HEALPix pixel */
    double x, y, mu, theta, phi;
    float value;
    for(int i = 0; i < w; i++) {
        for(int j = 0; j < h; j++) {
            value = 1. - image(i+xmargin,j+ymargin,0,0)/255.;

            /* Sample each pixel multiple times */
            for(int is = 0; is < nsamp; is++) {
                for(int js = 0; js < nsamp; js++) {
                    x = -1 + 2.*(i + (is+0.5)/nsamp)/w;
                    y =  1 - 2.*(j + (js+0.5)/nsamp)/h;
                    if(inverse_hammer(x, y, &mu, &phi) == 0) {
                        theta = acos(mu);
                        ang2pix_nest(nside, theta, phi, &ipix);
                        map[ipix] += value;
                        count[ipix] += 1;
                    }
                }
            }
        }
    }

    /* Take average map values over each HEALPix pixel */
    for(ipix = 0; ipix < npix; ipix++) {
        if(count[ipix] == 0) {
            /* HEALPix pixel didn't get updated, image pixels must be large */
            pix2ang_nest(nside, ipix, &theta, &phi);
            hammer(cos(theta), phi, &x, &y);
            map[ipix] = image(w*(1+x)/2. + xmargin, h*(1-y)/2. + ymargin, 0);
        }
        else
            map[ipix] /= count[ipix];
    }

    /* Write map to file, overwriting any previous file by the same name */
    remove(hpxfile);    // ignore errors
    write_healpix_map(map, nside, hpxfile, 1, "E");

    free(map);
    free(count);
    return 0;
}
