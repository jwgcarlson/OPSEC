/* multipole
 *
 * Computes monopole (l=0), quadrupole (l=2), and hexadecapole (l=4) moments
 * of the 2-point correlation function $\xi_\ell(r)$. */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>
using std::vector;

#include "Spline.h"

/* Spherical Bessel functions $j_n(x)$ for 0 <= n <= 4. */
static double sj0(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? 1 - x2*(1/6. + x2*(1/120. - x2*(1/5040.)))
                       : sin(x)/x;
}
static double sj1(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x*(1/3. - x2*(1/30. + x2*(1/840.)))
                       : -(x*cos(x) - sin(x))/x2;
}
static double sj2(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x2*(1/15. - x2*(1/210. + x2*(1/7560.)))
                       : -(3.*x*cos(x) + (x2 - 3.)*sin(x))/(x*x2);
}
static double sj3(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x*x2*(1/105. - x2*(1/1890.))
                       : (x*(x2 - 15.)*cos(x) - (6.*x2 - 15.)*sin(x))/(x2*x2);
}
static double sj4(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x2*x2*(1/945. - x2*(1/20790.))
                       : (x*(10.*x2 - 105.)*cos(x) + (x2*(x2 - 45.) + 105.)*sin(x))/(x*x2*x2);
}

/* Compute the quantity
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * using direct Simpson integration. */
static void ComputeXiLM(int l, int m, double kmin, double kmax, Spline P, int Nr, const double r[], double xi[]) {
    const int Nk = 65536;     // number of _subintervals_, not number of points; must be even
    int i, j, mult;
    double k, dk = (kmax - kmin)/Nk;

    /* Choose appropriate spherical Bessel function */
    double (*sj)(double x);
    if(l == 0)      sj = sj0;
    else if(l == 1) sj = sj1;
    else if(l == 2) sj = sj2;
    else if(l == 3) sj = sj3;
    else if(l == 4) sj = sj4;
    else {
        fprintf(stderr, "ComputeXiLM: l = %d not supported\n", l);
        return;
    }

    for(i = 0; i < Nr; i++)
        xi[i] = 0;

    /* Integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
    #pragma omp parallel for private(i,k,mult)
    for(j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 4 - 2*(j % 2) - 3*(j == 0) - 3*(j == Nk);
        mult *= P(k) * pow(k, m);
        for(i = 0; i < Nr; i++)
            #pragma omp atomic
            xi[i] += mult * sj(k*r[i]);
    }
    for(i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);
}

int main(int argc, char* argv[]) {
    const int Nr = 1024;
    const double rmin = 1;
    const double rmax = 200;
    const double kmin = 0.;
    const double kmax = 10.;

    if(argc != 5) {
        fprintf(stderr, "Usage: %s pkfile f b outfile\n", argv[0]);
        return 1;
    }

    /* Load power spectrum */
    const char* pkfile = argv[1];
    Spline P = LinearSpline(pkfile);
    if(!P) {
        fprintf(stderr, "Error reading P(k)n\n");
        return 1;
    }

    double f = strtod(argv[2], NULL);
    double b = strtod(argv[3], NULL);
    if(f <= 0 || b <= 0) {
        fprintf(stderr, "Invalid values, f = %g, b = %g\n", f, b);
        return 1;
    }

    const char* outfile = argv[4];
    FILE* fout = fopen(outfile, "w");
    if(!fout) {
        fprintf(stderr, "Could not open %s for writing\n", outfile);
        return 1;
    }

    /* Compute multipoles */
    vector<double> r(Nr), xi0(Nr), xi2(Nr), xi4(Nr);
    for(int i = 0; i < Nr; i++)
        r[i] = rmin + i*(rmax - rmin)/(Nr-1);
    ComputeXiLM(0, 2, kmin, kmax, P, Nr, &r[0], &xi0[0]);
    ComputeXiLM(2, 2, kmin, kmax, P, Nr, &r[0], &xi2[0]);
    ComputeXiLM(4, 2, kmin, kmax, P, Nr, &r[0], &xi4[0]);

    /* Multiply by appropriate Kaiser factors */
    double beta = f/b;
    double A0 = 1 + (2/3.)*beta + (1/5.)*beta*beta;
    double A2 = (4/3.)*beta + (4/7.)*beta*beta;
    double A4 = (8/35.)*beta*beta;
    for(int i = 0; i < Nr; i++) {
        xi0[i] *=  A0;
        xi2[i] *= -A2;
        xi4[i] *=  A4;
    }

    /* Write results to file */
    fprintf(fout, "# xi_ell(r) for f = %g, b = %g, with P(k) from %s\n", f, b, pkfile);
    for(int i = 0; i < Nr; i++)
        fprintf(fout, "%e %e %e %e\n", r[i], xi0[i], xi2[i], xi4[i]);
    fclose(fout);

    return 0;
}
