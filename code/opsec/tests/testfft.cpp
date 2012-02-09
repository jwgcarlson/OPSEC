/* testfft
 *
 * Test program for C++ implementation of FFTLog algorithm.  Returns 0 on
 * success, 1 on failure. */

#include <cmath>
#include <cstdio>

#include "CuTest.h"
#include "fftlog.h"

static inline double sq(double x) { return x*x; }

/* Compute the \xi(r) corresponding to a Gaussian,
 *   P(k) = (2\pi)^{3/2} e^{-k^2/}.
 * The answer should be \xi(r) = e^{-r^2/2} to reasonable precision. */
void TestGaussian(CuTest* tc) {
    /* Define P(k) */
    const int N = 8192;
    const double kmin = 1e-5, kmax = 1e5;
    double k[N], pk[N];
    for(int i = 0; i < N; i++) {
        k[i] = kmin * exp(i*log(kmax/kmin)/(N-1));
        pk[i] = pow(2*M_PI, 1.5) * exp(-k[i]*k[i]/2);
    }

    /* Compute \xi(r) */
    double r[N], xi[N];
    pk2xi(N, k, pk, r, xi);
//    printf("rmin = %g, rmax = %g\n", r[0], r[N-1]);

    /* Compute the integral of the square of the function and its residual */
    double integral = 0;
    double residual = 0;
    for(int i = 0; i < N-1; i++) {
        integral += (r[i+1] - r[i]) * 0.5 * ( exp(-r[i+1]*r[i+1]) + exp(-r[i]*r[i]) );
        residual += (r[i+1] - r[i]) * 0.5 * ( sq(xi[i+1] - exp(-r[i+1]*r[i+1]/2)) + sq(xi[i] - exp(-r[i]*r[i]/2)) );
    }

//    for(int i = 0; i < N; i++)
//        printf("%g %g %g\n", r[i], xi[i], exp(-r[i]*r[i]/2));

    /* The function itself should integrate to \sqrt{\pi}/2 */
//    printf("integral = %g, sqrt(pi)/2 = %g\n", integral, sqrt(M_PI)/2);
    CuAssertTrue(tc, fabs(integral - sqrt(M_PI)/2) < 1e-3);

    /* The residual should be within 0.1% */
//    printf("residual = %g\n", residual);
    CuAssertTrue(tc, residual/integral < 1e-3);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestGaussian);

    return suite;
}

int main(int argc, char* argv[]) {
    CuString* output = CuStringNew();
    CuSuite* suite = CuSuiteNew();

    CuSuiteAddSuite(suite, GetSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);

    return (suite->failCount == 0) ? 0 : 1;
}
