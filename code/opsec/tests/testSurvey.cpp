/* testSurvey
 *
 * Test Survey base class. */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Cell.h"
#include "Survey.h"
#include "cfg.h"

#include "CuTest.h"

void TestSurveyConstructor(CuTest* tc) {
    Config cfg = cfg_new();
    cfg_set(cfg, "XMin", "0");
    cfg_set(cfg, "XMax", "1");
    cfg_set(cfg, "YMin", "2");
    cfg_set(cfg, "YMax", "3");
    cfg_set(cfg, "ZMin", "4");
    cfg_set(cfg, "ZMax", "5");
    cfg_set(cfg, "Nx", "6");
    cfg_set(cfg, "Ny", "7");
    cfg_set(cfg, "Nz", "8");
    Survey* survey = new Survey(CoordSysCartesian, cfg);

    CuAssertTrue(tc, survey->GetCoordinateSystem() == CoordSysCartesian);
    CuAssertTrue(tc, survey->XMin == 0);
    CuAssertTrue(tc, survey->XMax == 1);
    CuAssertTrue(tc, survey->YMin == 2);
    CuAssertTrue(tc, survey->YMax == 3);
    CuAssertTrue(tc, survey->ZMin == 4);
    CuAssertTrue(tc, survey->ZMax == 5);
    CuAssertTrue(tc, survey->Nx == 6);
    CuAssertTrue(tc, survey->Ny == 7);
    CuAssertTrue(tc, survey->Nz == 8);

    delete survey;
    cfg_destroy(cfg);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestSurveyConstructor);

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
