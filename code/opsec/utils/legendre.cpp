/* legendre
 *
 * Compute the Legendre coefficients of the 2-point correlation function xi(r)
 * evaluated on a thin spherical shell.
 *
 *   f(x) = \xi(R\sqrt{2(1-x)})
 *        = \sum_n a_n P_n(x)
 *    a_n = \frac{2n+1}{2} \int_{-1}^1 f(x) P_n(x) dx
 *
 * For our purposes we want to know $a_n/(2n+1)$. */

#include <cmath>
#include <cstdio>
#include <vector>

#include "Model.h"
#include "QSpline.h"

/* Number of Legendre coefficients to compute */
int numcoeffs = 1024;

/* Number of grid points to use on the interval [-1,1] */
int N = 8192;

/* Radius of shell */
double R = 500.;

const char* pkfile = "pk.dat";

int n1 = 512;
int n2 = 256;

double r1 = 200.;
double r2 = 2000.;

int main(int argc, char* argv[]) {
    /* Parse command line arguments */
    // TODO

    /* Load P(k) */
    Spline P = CubicSpline(pkfile);

    /* Comptue xi(r) */
    std::vector<double> ar = QSpline::MakeArray(n1, n2, r1, r2);
    std::vector<double> axi(n1 + n2);
    ComputeXiLM(0, 2, P, n1 + n2, &ar[0], &axi[0]);
    QSpline xi(n1, n2, ar, axi);

    /* Evaluate function f(x) on a grid */
    std::vector<double> ax(N+1);
    std::vector<double> af(N+1);
    for(int i = 0; i <= N; i++) {
        ax[i] = -1. + i*2./N;
        af[i] = xi(R*sqrt(2*(1 - ax[i])));
    }

    /* Array of result values */
    std::vector<double> a(numcoeffs, 0);

    double dx = 2./N;
    for(int i = 0; i <= N; i++) {
        double mult = (2 + 2*(i % 2) - (i == 0) - (i == N)) * dx/3.;
        double x = -1. + i*2./N;
        double f = xi(R*sqrt(2*(1 - x)));

        double p0 = 1.;
        a[0] += 0.5 * mult * f * p0;

        double p1 = x;
        a[1] += 1.5 * mult * f * p1;

        /* Upward recursion for Legendre polynomials */
        for(int n = 2; n < numcoeffs; n++) {
            double pn = ( (2*n-1)*x*p1 - (n-1)*p0 ) / n;
            a[n] += (n + 0.5) * mult * f * pn;
            p0 = p1;
            p1 = pn;
        }
    }

    /* Print results */
    printf("#  n   a_n   a_n/(2n+1)\n");
    for(int n = 0; n < numcoeffs; n++)
        printf("% 4d %8e %8e\n", n, a[n], a[n]/(2*n+1));

    return 0;
}
