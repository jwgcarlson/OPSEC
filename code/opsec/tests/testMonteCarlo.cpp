/* testMonteCarlo
 *
 * Test program for MonteCarloIntegral class. Uses the CuTest unit testing
 * framework.  Returns 0 on success, 1 on failure. */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "CuTest.h"
#include "MonteCarlo.h"
#include "rng.h"


class Integral1 : public MonteCarloIntegral {
public:
    /* (1-cos(1)) * sin(1) * (exp(1)-1) = 0.664669679978137714 */
    static const double exact = 0.664669679978137714;

    Integral1() : MonteCarloIntegral(3) {
        maxeval = 100000;
    }

    void Integrand(const double x[], double* f, double*) {
        f[0] = sin(x[0]) * cos(x[1]) * exp(x[2]);
    }
};

void TestIntegral1(CuTest* tc) {
    int neval;
    double value, error;
    Integral1 I1;
    I1.Integrate(&value, NULL, &error, &neval, 0, 0.001);
    printf("I1 = %.12f +/- %.12f (%d evaluations, exact value %.12f)\n", value, error, neval, I1.exact);
    CuAssertDblEquals(tc, I1.exact, value, 2*error);

    int maxeval = I1.maxeval;
    Sobol sobol;
    SobolIni(&sobol, maxeval, 3);
    double* points = (double*) malloc(maxeval*3*sizeof(double));
    for(int k = 0; k < maxeval; k++) {
        SobolGet(&sobol, &points[3*k]);
//        printf("x[%d] = %g %g %g\n", k, points[3*k], points[3*k+1], points[3*k+2]);
    }
    I1.Integrate(&value, points, NULL, &error, &neval, 0, 0.001);
    printf("I1 = %.12f +/- %.12f (%d evaluations, exact value %.12f)\n", value, error, neval, I1.exact);
    CuAssertDblEquals(tc, I1.exact, value, 2*error);

    srand(1234);
    for(int i = 0; i < 3*maxeval; i++) {
        points[i] = rand()/((double)RAND_MAX+1);
    }
    I1.Integrate(&value, points, NULL, &error, &neval, 0, 0.001);
    printf("I1 = %.12f +/- %.12f (%d evaluations, exact value %.12f)\n", value, error, neval, I1.exact);
    CuAssertDblEquals(tc, I1.exact, value, 2*error);

    free(points);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestIntegral1);

    return suite;
}

int main(int argc, char **argv) {
    CuString* output = CuStringNew();
    CuSuite* suite = CuSuiteNew();

    CuSuiteAddSuite(suite, GetSuite());

    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);

    return (suite->failCount == 0) ? 0 : 1;
}
