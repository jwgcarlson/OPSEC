/* realmock
 *
 * Reads in particle data from an N-body simulation or a Gaussian mock
 * generator, and produces a more realistic mock catalog according to a given
 * survey specification. */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <vector>
using std::vector;

#include "Spline.h"
#include "Survey.h"
#include "abn.h"
#include "cfg.h"
#include "chealpix.h"
#include "rng.h"

const char* usage =
    "Usage: realmock blah balh blalh (TODO)\n";

const char* required_cfgopts =
    "input_galfile,input_coordsys,output_galfile,output_coordsys,maskfile,radialfile,output_ngals,L,RMax";

/* Globals */
int output_ngals = 0;
const char* maskfile = NULL;
const char* radialfile = NULL;
float* mask = NULL;
Spline radial;
double L = 0.;
double origin[3] = { 0., 0., 0. };
double RMin = -1;
double RMax = -1;
const char* input_galfile = NULL;
const char* output_galfile = NULL;
const char* input_coordsys = NULL;
const char* output_coordsys = NULL;
unsigned int seed = 1337;

enum CoordSys {
    SPHERICAL = 0,
    CARTESIAN = 1,
    UNDEFINED
};

float read_mask(const float* mask, long nside, int nest, double theta, double phi) {
    long ipix;
    if(nest)
        ang2pix_nest(nside, theta, phi, &ipix);
    else
        ang2pix_ring(nside, theta, phi, &ipix);
    return mask[ipix];
}

int main(int argc, char* argv[]) {
    Config cfg = cfg_new();
    Config input_galopts = cfg_new();
    Config output_galopts = cfg_new();

    /* Parse command line switches */
    int opt;
    while((opt = getopt(argc, argv, "hc:")) != -1) {
        switch(opt) {
        case 'h':
            fputs(usage, stdout);
            exit(0);
            break;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fputs(usage, stderr);
            exit(1);
            break;
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Debugging... */
    printf("# Config options\n");
    cfg_write(cfg, stdout);
    printf("\n");

    /* Read configuration options */
    if(!cfg_has_keys(cfg, required_cfgopts, ",")) {
        fprintf(stderr, "realmock: missing required configuration options\n");
        exit(1);
    }
    output_ngals = cfg_get_int(cfg, "output_ngals");
    maskfile = cfg_get(cfg, "maskfile");
    radialfile = cfg_get(cfg, "radialfile");
    L = cfg_get_double(cfg, "L");
    if(cfg_has_key(cfg, "origin"))
        cfg_get_array_double(cfg, "origin", 3, origin);
    if(cfg_has_key(cfg, "RMin"))
        RMin = cfg_get_double(cfg, "RMin");
    if(cfg_has_key(cfg, "RMax"))
        RMax = cfg_get_double(cfg, "RMax");
    input_galfile = cfg_get(cfg, "input_galfile");
    output_galfile = cfg_get(cfg, "output_galfile");
    input_coordsys = cfg_get(cfg, "input_coordsys");
    output_coordsys = cfg_get(cfg, "output_coordsys");
    if(cfg_has_key(cfg, "seed"))
        seed = cfg_get_uint(cfg, "seed");

    /* Read HEALPix mask */
    long nside;
    char coordsys[9], ordering[32];
    int nest;
    mask = read_healpix_map(maskfile, &nside, coordsys, ordering);
    if(strncmp(ordering, "NESTED", 6) == 0)
        nest = 1;
    else if(strncmp(ordering, "RING", 4) == 0)
        nest = 0;
    else {
        fprintf(stderr, "realmock: unrecognized HEALPix ordering: %s\n", ordering);
        exit(1);
    }

    /* Read radial selection function */
    radial = LinearSpline(radialfile);

    /* Guess RMin and RMax if they are not specified */
    if(RMin < 0 || RMax < 0) {
        double r = 0.;
        double dr = 10.;        // 10 Mpc/h steps
        while(radial(r) <= 0) {
            r += dr;
            if(r > 100000) {    // something is wrong
                fprintf(stderr, "realmock: RMin/RMax not specified and radial number density appears to be zero\n");
                exit(1);
            }
        }
        if(RMin < 0)
            RMin = r - dr;
        if(RMax < 0) {
            while(radial(r) > 0) // this assumes the support of the radial number density is connected
                r += dr;
            RMax = r;
        }
    }

    /* Read input particle file header */
    FILE* input_fgal = fopen(input_galfile, "r");
    if(input_fgal == NULL) {
        fprintf(stderr, "realmock: could not read galaxies from %s\n", input_galfile);
        exit(1);
    }
    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    if(abn_read_header(input_fgal, &n, &size, &endian, fmt, input_galopts) != 0
       || size != 4*sizeof(float)
       || strcmp(fmt, "4f") != 0)
    {
        fprintf(stderr, "realmock: error reading particle file header from %s\n", input_galfile);
        exit(1);
    }

    int input_ngals = (int) n;
    double input_nbar = input_ngals / (L*L*L);

    /* Compute normalizing selection factor to get input_ngals ~= output_ngals.
     * That is, treat the spline 'radial' as giving only the _shape_ of the
     * radial selection function, and compute the appropriate factor Q so that
     * Q*radial(r) gives the actual radial number density $\bar{n}(r)$. */
    double Q;
    {
        /* Compute radial integral $\int \bar{n}(r) r^2 dr$ using trapezoid rule */
        double r, dr = (RMax - RMin)/8192;
        double radial_integral = 0.5*RMin*RMin*radial(RMin) + 0.5*RMax*RMax*radial(RMax);
        for(int i = 1; i < 8192; i++) {
            r = RMin + i*dr;
            radial_integral += r*r*radial(r);
        }
        radial_integral *= dr;
        if(radial_integral <= 0.) {
            fprintf(stderr, "realmock: radial integral = %g\n", radial_integral);
            exit(1);
        }

        double angular_integral = 0.;
        long npix = 12*nside*nside;
        for(long ipix = 0; ipix < npix; ipix++)
            angular_integral += mask[ipix];
        angular_integral *= 4*M_PI/npix;
        if(radial_integral <= 0.) {
            fprintf(stderr, "realmock: angular integral = %g\n", angular_integral);
            exit(1);
        }

        Q = output_ngals / (radial_integral * angular_integral);
        printf("Q = %g\n", Q);
    }

    /* Compare maximum requested number density to available number density */
    double max_requested_nbar = 0.;
    for(int i = 0; i <= 8192; i++) {
        double nbar = Q * radial(i*RMax/8192.);
        if(nbar > max_requested_nbar)
            max_requested_nbar = nbar;
    }
    if(max_requested_nbar > input_nbar) {
        fprintf(stderr, "realmock: can't produce %d output galaxies with given radial profile\n", output_ngals);
        exit(1);
    }

    /* Open output file to make sure we can write to it */
    FILE* output_fgal = fopen(output_galfile, "w");
    if(output_fgal == NULL) {
        fprintf(stderr, "realmock: cannot write to %s\n", output_galfile);
        exit(1);
    }

    /* Initialize random number generator */
    rng_init(seed);

    /* Decide on coordinate systems for input particles and output galaxies */
    int input_coordsys_flag, output_coordsys_flag;
    if(strcmp(input_coordsys, "file") == 0)
        input_coordsys = cfg_get(input_galopts, "coordsys");
    if(strcmp(input_coordsys, "spherical") == 0)
        input_coordsys_flag = SPHERICAL;
    else if(strcmp(input_coordsys, "cartesian") == 0)
        input_coordsys_flag = CARTESIAN;
    else {
        fprintf(stderr, "realmock: invalid or missing input_coordsys option: %s\n", input_coordsys);
        exit(1);
    }
    if(strcmp(output_coordsys, "spherical") == 0)
        output_coordsys_flag = SPHERICAL;
    else if(strcmp(output_coordsys, "cartesian") == 0)
        output_coordsys_flag = CARTESIAN;
    else {
        fprintf(stderr, "realmock: invalid or missing output_coordsys option: %s\n", output_coordsys);
        exit(1);
    }

    /* Start sampling galaxies */
    vector<Galaxy> output_galaxies;
    Galaxy g;
    double r, mu, phi;
    double x, y, z;
    bool radial_pass, angular_pass;
    int nread;
    for(int i = 0; i < input_ngals; i++) {
        nread = fread(&g, sizeof(Galaxy), 1, input_fgal);
        if(nread != 1) {
            fprintf(stderr, "realmock: error reading particle number %d: nread = %d\n", i, nread);
            exit(1);
        }

        /* Get spherical coordinates for particle */
        switch(input_coordsys_flag) {
        case SPHERICAL:
            r = g.r;
            mu = g.mu;
            phi = g.phi;
            break;
        case CARTESIAN:
            x = g.x - origin[0];
            y = g.y - origin[1];
            z = g.z - origin[2];
            r = sqrt(x*x + y*y + z*z);
            mu = z/r;
            phi = M_PI + atan2(y, x);
            break;
        }

        /* Decide whether or not to keep it as a galaxy */
        radial_pass = r >= RMin &&
                      r <= RMax &&
                      input_nbar * rng_uniform() < Q * radial(r);
        angular_pass = rng_uniform() < read_mask(mask, nside, nest, acos(mu), phi);

        if(radial_pass && angular_pass) {
            switch(output_coordsys_flag) {
            case SPHERICAL:
                g.r = r;
                g.mu = mu;
                g.phi = phi;
                break;
            case CARTESIAN:
                g.x = r*sqrt(1-mu*mu)*cos(phi);
                g.y = r*sqrt(1-mu*mu)*sin(phi);
                g.z = r*mu;
                break;
            }
            output_galaxies.push_back(g);
        }
    }

    /* Write output galaxies to file */
    fprintf(output_fgal, "# Mock galaxy catalog created by realmock\n");
    cfg_set(output_galopts, "coordsys", output_coordsys);
    cfg_set(output_galopts, "input_galfile", input_galfile);
    cfg_set(output_galopts, "maskfile", maskfile);
    cfg_set(output_galopts, "radialfile", radialfile);
    cfg_set_uint(output_galopts, "seed", seed);
    abn_write(output_fgal, &output_galaxies[0], output_galaxies.size(), "4f", output_galopts);
    fclose(output_fgal);

    cfg_destroy(input_galopts);
    cfg_destroy(output_galopts);
    free(mask);
    return 0;
}
