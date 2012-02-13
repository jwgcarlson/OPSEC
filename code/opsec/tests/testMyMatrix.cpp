/* testMyMatrix
 *
 * Test program for MyMatrix and MyVector wrapper classes. */

#include <cmath>
#include <cstdlib>

#include "CuTest.h"
#include "MyMatrix.hpp"


void TestMyVector(CuTest* tc) {
    const int n = 128;
    double buf[n];

    MyVector<double> x(n, 1, buf);
    for(int i = 0; i < n; i++)
        x[i] = sqrt(i);
    CuAssertDblEquals(tc, x[10], sqrt(10), 1e-10);
    CuAssertDblEquals(tc, x.norm(), sqrt(n*(n-1)/2.), 1e-10);

    /* Vectors just point to the underlying memory, so y and x both refer to buf. */
    MyVector<double> y = x;
    y[10] = 300;
    CuAssertDblEquals(tc, x[10], 300, 1e-10);
    y[10] = sqrt(10);

    /* A copy, however, points to a new memory region. */
    double buf2[n];
    MyVector<double> z = x.copy(buf2);
    CuAssertDblEquals(tc, z[10], sqrt(10), 1e-10);
    CuAssertDblEquals(tc, z.norm(), sqrt(n*(n-1)/2.), 1e-10);
    z[10] = 300;
    CuAssertDblEquals(tc, x[10], sqrt(10), 1e-10);
}

void TestMyMatrix(CuTest* tc) {
    const int m = 32;
    const int n = 32;
    double buf[m*n];

    MyMatrix<double> A(m, n, ColumnMajor, 0, buf);
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
            A(i,j) = i+j;
    CuAssertDblEquals(tc, A(10,10), 20, 1e-10);
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestMyVector);
    SUITE_ADD_TEST(suite, TestMyMatrix);

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
