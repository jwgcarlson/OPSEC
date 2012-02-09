
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <unistd.h>

#include "abn.h"
#include "cfg.h"
#include "Survey.h"

template<typename T>
static inline T cb(T x) { return x*x*x; }

/* Generate a random number uniformly in the half-open interval [0,1). */
static double uniform() {
    return rand() / (RAND_MAX + 1.);
}

int main(int argc, char* argv[]) {
    Config cfg = cfg_new();

    /* Parse command line switches */
    int opt;
    const char* optstring = "hc:";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            printf("Usage: %s [-h] [-c <config>] [CONFIG OPTIONS]\n", argv[0]);
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fprintf(stderr, "Usage: %s [-h] [-c <config>] [CONFIG OPTIONS]\n", argv[0]);
            return 1;
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    if(!cfg_has_keys(cfg, "galfile,Ngals,RMin,RMax,MuMin,MuMax,PhiMin,PhiMax", ",")) {
        fprintf(stderr, "Must provide config options galfile, Ngals, and {R,Mu,Phi}{Min,Max}\n");
        return 1;
    }

    const char* galfile = cfg_get(cfg, "galfile");
    int Ngals = cfg_get_int(cfg, "Ngals");
    double RMin = cfg_get_double(cfg, "RMin");
    double RMax = cfg_get_double(cfg, "RMax");
    double MuMin = cfg_get_double(cfg, "MuMin");
    double MuMax = cfg_get_double(cfg, "MuMax");
    double PhiMin = cfg_get_double(cfg, "PhiMin");
    double PhiMax = cfg_get_double(cfg, "PhiMax");

    if(cfg_has_key(cfg, "seed"))
        srand(cfg_get_uint(cfg, "seed"));

    double V = (cb(RMax) - cb(RMin))/3 * (MuMax - MuMin) * (PhiMax - PhiMin);
    printf("density Ngals/V = %g\n", Ngals/V);

    /* Randomly sprinkle galaxies within the specified region. */
    Galaxy* gals = (Galaxy*) malloc(Ngals*sizeof(Galaxy));
    for(int g = 0; g < Ngals; g++) {
        gals[g].r = cbrt(cb(RMin) + uniform()*(cb(RMax) - cb(RMin)));
        gals[g].mu = MuMin + uniform()*(MuMax - MuMin);
        gals[g].phi = PhiMin + uniform()*(PhiMax - PhiMin);
        gals[g].w = 1;
    }

    /* Write galaxies to file */
    FILE* fgals = fopen(galfile, "w");
    abn_write(fgals, (void*)gals, (size_t)Ngals, "3f", NULL);
    fclose(fgals);

    free(gals);
    return 0;
}
