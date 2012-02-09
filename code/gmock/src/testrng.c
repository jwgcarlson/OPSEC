
#include <math.h>
#include <stdio.h>

#include "rng.h"

double fact(int k) {
    double f = 1.;
    while(k > 0)
        f *= (k--);
    return f;
}

int main() {
    int i;
    unsigned int seed = 1776;

    rng_init(seed);

    /* Test uniform random numbers */
    {
        FILE* funif = fopen("uniform.dat", "w");
        const int npoints = 1000000;
        const int nbins = 128;
        double bins[nbins];
        int i, b;
        double x;
        for(b = 0; b < nbins; b++)
            bins[b] = 0;
        for(i = 0; i < npoints; i++) {
            x = rng_uniform();
            bins[(int)(x*nbins)] += 1;
        }

        for(b = 0; b < nbins; b++)
            fprintf(funif, "%5.5f %10.10f\n", (b+0.5)/nbins, nbins*bins[b]/npoints);
        fclose(funif);
    }

    /* Test normally distributed random numbers */
    {
        FILE* fnorm = fopen("normal.dat", "w");
        const int npoints = 1000000;
        const int nbins = 128;
        double bins[nbins];
        int i, b;
        double x;
        for(b = 0; b < nbins; b++)
            bins[b] = 0;
        for(i = 0; i < npoints; i++) {
            x = rng_normal();
            if(x >= -4 && x < 4)
                bins[(int)(nbins*(x+4)/8)] += 1;
        }

        for(b = 0; b < nbins; b++)
            fprintf(fnorm, "%5.5f %10.10f\n", -4 + 8*(b+0.5)/nbins, nbins*bins[b]/(8*npoints));
        fclose(fnorm);
    }

    /* Test Poisson sampling function */
    {
        FILE* fpois = fopen("poisson.dat", "w");
        const int npoints = 1000000;
        const int kmax = 128;
        const double lambda = 50.;
        double c[kmax];
        int i, k;
        for(k = 0; k < kmax; k++)
            c[k] = 0;
        for(i = 0; i < npoints; i++) {
            k = (int)rng_poisson(lambda);
            if(k < kmax) c[k] += 1;
        }

        for(k = 0; k < kmax; k++)
            fprintf(fpois, "%d %10.10f %10.10f\n", k, c[k]/npoints, exp(-lambda)*pow(lambda, k)/fact(k));
        fclose(fpois);
    }

#if 0
    /* Track down source of nan's in gcorr */
    {
        int N = 1024;
        int ix, iy, iz;
        double u1, u2;
        rng_init(1);
        for(ix = 0; ix < N; ix++) {
            for(iy = 0; iy < N; iy++) {
                for(iz = 0; iz <= N/2; iz++) {
                    u1 = rng_uniform();
                    u2 = rng_uniform();

                    /* u1 = 1 for both cases !!! */
                    if(ix == 313 && iy == 722 && iz == 210)
                        printf("(313,722,210): u1 = %g, u2 = %g\n", u1, u2);
                    if(ix == 848 && iy == 940 && iz == 350)
                        printf("(848,940,350): u1 = %g, u2 = %g\n", u1, u2);

                    if(u1 == 0 || u2 == 0)
                        printf("zero!\n");
                }
            }
        }
    }
#endif

    return 0;
}
