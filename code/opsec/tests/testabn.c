/* testabn
 *
 * Test program for abn_* routines.  Uses the CuTest unit testing framework.
 * Returns 0 on success, 1 on failure. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CuTest.h"
#include "abn.h"
#include "cfg.h"


void TestAbnReadWrite(CuTest* tc) {
    const size_t ndata = 100;
    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    double *data;
    int i;
    FILE* fp;

    fp = tmpfile();

    data = (double*)malloc(ndata*sizeof(double));
    for(i = 0; i < ndata; i++)
        data[i] = cos(i);
    abn_write(fp, data, ndata, "d", NULL);

    free(data);
    data = NULL;
    rewind(fp);

    abn_read(fp, (void**)&data, &n, &size, &endian, fmt, NULL);
    CuAssertTrue(tc, n == ndata);
    CuAssertTrue(tc, size == sizeof(double));
    CuAssertTrue(tc, endian == abn_endian());
    CuAssertStrEquals(tc, fmt, "d");
    for(i = 0; i < ndata; i++)
        CuAssertTrue(tc, fabs(data[i] - cos(i)) < 1e-15);

    free(data);
    fclose(fp);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestAbnReadWrite);

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
