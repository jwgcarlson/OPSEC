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

    MyVector<double> x(n, buf);
    for(int i = 0; i < n; i++)
        x[i] = sqrt(i);
    CuAssertDblEquals(tc, sqrt(10), x[10], 1e-10);
    CuAssertDblEquals(tc, sqrt(n*(n-1)/2.), x.norm(), 1e-10);

    /* Vectors just point to the underlying memory, so y and x both refer to buf. */
    MyVector<double> y = x;
    y[10] = 300;
    CuAssertDblEquals(tc, 300, x[10], 1e-10);
    y[10] = sqrt(10);

    /* A copy, however, points to a new memory region. */
    double buf2[n];
    MyVector<double> z = x.newcopy(buf2);
    CuAssertDblEquals(tc, sqrt(10), z[10], 1e-10);
    CuAssertDblEquals(tc, sqrt(n*(n-1)/2.), z.norm(), 1e-10);
    z[10] = 300;
    CuAssertDblEquals(tc, sqrt(10), x[10], 1e-10);
}

void TestMyMatrix(CuTest* tc) {
    const int m = 32;
    const int n = 32;
    double bufa[m*n];

    MyMatrix<double> A(m, n, bufa);
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
            A(i,j) = i - j;
    CuAssertDblEquals(tc, 10, A(15,5), 1e-10);

    MyMatrix<double> At = A.transpose();
    CuAssertDblEquals(tc, -10, At(15,5), 1e-10);

    double bufb[3*3] = { 1, 2, 3,
                         4, 5, 6,
                         7, 8, 9 };
    double bufx[3] = { 1, 2, 3 };
    double bufy[3];
    A = MyMatrix<double>(3, 3, bufb);
    MyVector<double> x(3, bufx);
    MyVector<double> y(3, bufy);
    multiply(A, x, y);
    CuAssertDblEquals(tc, 14, y[0], 1e-10);
    CuAssertDblEquals(tc, 32, y[1], 1e-10);
    CuAssertDblEquals(tc, 50, y[2], 1e-10);

    CuAssertDblEquals(tc, 15, A.trace(), 1e-10);
}

void TestMatrixSqrt(CuTest* tc) {
    /* Explicit example borrowed from Matlab documentation */
    double xvalues[5*5] = { 5, -4,  1,  0,  0,
                           -4,  6, -4,  1,  0,
                            1, -4,  6, -4,  1,
                            0,  1, -4,  6, -4,
                            0,  0,  1, -4,  5 };
    double yvalues[5*5];

    MyMatrix<double> X(5, 5, xvalues);
    MyMatrix<double> Y(5, 5, yvalues);
    sqrtm(X, Y);

#define CK(i, j, ex) CuAssertDblEquals(tc, ex, Y(i-1,j-1), 1e-10);
    CK(1,1, 2.)  CK(1,2, -1.) CK(1,3, 0.)  CK(1,4, 0.)  CK(1,5, 0.)
    CK(2,1, -1.) CK(2,2, 2.)  CK(2,3, -1.) CK(2,4, 0.)  CK(2,5, 0.)
    CK(3,1, 0.)  CK(3,2, -1.) CK(3,3, 2.)  CK(3,4, -1.) CK(3,5, 0.)
    CK(4,1, 0.)  CK(4,2, 0.)  CK(4,3, -1.) CK(4,4, 2.)  CK(4,5, -1.)
    CK(5,1, 0.)  CK(5,2, 0.)  CK(5,3, 0.)  CK(5,4, -1.) CK(5,5, 2.)
#undef CK
}

void TestMatrixInverse(CuTest* tc) {
    double avalues[3*3] = { 1,  1,  2,
                            2,  4, -3,
                            3,  6, -5 };
    double bvalues[3*3];

    MyMatrix<double> A(3, 3, avalues);
    MyMatrix<double> B(3, 3, bvalues);
    inv(A, B);

#define CK(i, j, ex) CuAssertDblEquals(tc, ex, B(i-1,j-1), 1e-10);
    CK(1,1, 2.)  CK(1,2, -17.) CK(1,3, 11.)
    CK(2,1, -1.) CK(2,2, 11.)  CK(2,3, -7.)
    CK(3,1, 0.)  CK(3,2, 3.)   CK(3,3, -2.)
#undef CK
}

void TestPackedMatrixIndex(CuTest* tc) {
    int n = 5;
    PackedMatrixIndex gr = 0;
    int i = 0;
    int j = 0;
    for(int r = 0; r < n*(n+1)/2; ++r) {
//        printf("ap[%d] = A(%d,%d)\n", r, gr.i, gr.j);
        CuAssertIntEquals(tc, i, gr.i);
        CuAssertIntEquals(tc, j, gr.j);
        ++j;
        if(j > i) {
            ++i;
            j = 0;
        }
        ++gr;
    }
}

CuSuite* GetSuite() {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, TestMyVector);
    SUITE_ADD_TEST(suite, TestMyMatrix);
    SUITE_ADD_TEST(suite, TestMatrixSqrt);
    SUITE_ADD_TEST(suite, TestMatrixInverse);
    SUITE_ADD_TEST(suite, TestPackedMatrixIndex);

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
