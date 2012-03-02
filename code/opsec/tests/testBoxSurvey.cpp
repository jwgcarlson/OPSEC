/* testBoxSurvey
 *
 * Test BoxSurvey class. */

#include <cmath>
#include <cstdio>

#include "BoxSurvey.h"
#include "SelectionFunc.h"
#include "SeparationFunc.h"

#include "CuTest.h"


void TestSeparationFunc(CuTest* tc) {
    /* 1 million galaxies over [0,100]^3 */
    Survey* survey = new BoxSurvey("none", 1., 100.);

    CuAssertTrue(tc, survey->GetCoordinateSystem() == CoordSysCartesian);

    SeparationFunc sep = survey->GetSeparationFunction();
    Point p1 = { 10, 20, 30 };
    Point p2 = { 70, 80, 90 };
    CuAssertDblEquals(tc, sqrt(3*40*40), sep.r(p1, p2), 1e-10);
    p2.z = 50;
    CuAssertDblEquals(tc, sqrt(2*40*40 + 20*20), sep.r(p1, p2), 1e-10);

    delete survey;
}

void TestSelectionFunc(CuTest* tc) {
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestSeparationFunc);
    SUITE_ADD_TEST(suite, TestSelectionFunc);

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
