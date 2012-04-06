#ifndef QSPLINE_H
#define QSPLINE_H

#include <vector>

#include "Spline.h"

/* Combination of linearly-spaced spline F1(x) and logarithmically-spaced
 * spline F2(x).
 *
 * Returns:
 *  x < 0          0
 *  0 <= x < x1    F1(x)
 *  x1 <= x < x2   F2(x)
 *  x >= x2        F2(x2)
 */

class QSpline {
public:
    QSpline();
    QSpline(int n1, int n2, const double x[], const double y[]);
    QSpline(int n1, int n2, const std::vector<double>& x, const std::vector<double>& y);
    ~QSpline();

    double operator()(double x) const;

    /* Creates an array of length n1+n2, with the first n1 elements linearly
     * spaced over [0,x1], and the last n2 elements logarithmically spaced over
     * [x1,x2].  Element number 0 is equal to 0, element number n1 is equal to
     * x1, and the last element (number n1+n2-1) is equal to x2. */
    static std::vector<double> MakeArray(int n1, int n2, double x1, double x2);

protected:
    Spline F1;
    Spline F2;
    double x1;          // domain for linearly-spaced spline, [0,x1)
    double x2;          // domain for logarithmically-spaced spline, [x1,x2)
    double y2;          // F(x2)
};

#endif // QSPLINE_H
