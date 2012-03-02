/* testsig
 *
 * Test signal matrix computation (sig.cpp). */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "BoxSurvey.h"
#include "Cell.h"
#include "SeparationFunc.h"
#include "TestModel.h"
#include "XiFunc.h";
#include "cfg.h"
#include "sig.h"

#include "CuTest.h"

void TestConstantXi(CuTest* tc) {
    /* xi = 1 = const */
    Model* model = new TestModel(1, 1., 0.);
    XiFunc xi = model->GetXi();

    Survey* survey = new BoxSurvey("none", 1.0, 1000);
    SeparationFunc sep = survey->GetSeparationFunction();

    /* Create cubical cells */
    Cell c1 = { 0, 0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1., 1. };
    Cell c2 = { 1, 1, 1.0, 1.5, 1.0, 1.5, 1.0, 1.5, 1., 0.125 };

    int neval;
    double Q = ComputeSignalC(c1, c2, xi, sep, 1e-5, 1e-10, &neval);
    double S = c1.Nbar * c2.Nbar * Q;
    printf("TestConstantXi: neval = %d\n", neval);
    CuAssertDblEquals(tc, 0.125, S, 1e-12);

    delete model;
    delete survey;
}

void TestGaussianXi(CuTest* tc) {
    const double nbar = 1.;
    const double L = 10.;
    const double sigma = 1/M_SQRT2;
    const double epsrel = 1e-5;

    /* \xi = (2\pi\sigma^2)^{-3/2} e^{-r^2/2\sigma^2} */
    Model* model = new TestModel(1, pow(2*M_PI*sigma*sigma, -1.5), 1./(2*sigma*sigma));
    XiFunc xi = model->GetXi();

    Survey* survey = new BoxSurvey("none", nbar, L);
    SeparationFunc sep = survey->GetSeparationFunction();

    /* Create cubical cell */
    Cell c = { 0, 0, 0., L, 0., L, 0., L, L*L*L, nbar*L*L*L };

    int neval;
    double Q = ComputeSignalC(c, c, xi, sep, epsrel, 1e-12, &neval);
    printf("TestGaussianXi: neval = %d\n", neval);
    double S = c.Nbar * c.Nbar * Q;
    double Sexpected = nbar*nbar * pow(L*erf(L/(2*M_SQRT2*sigma)), 3);
    CuAssertDblEquals(tc, Sexpected, S, 5*epsrel * S);

    delete model;
    delete survey;
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestConstantXi);
    SUITE_ADD_TEST(suite, TestGaussianXi);

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
