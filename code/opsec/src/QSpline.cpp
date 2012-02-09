#include <cassert>
#include <cmath>

#include "QSpline.h"

vector<double> QSpline::MakeArray(int n1, int n2, double x1, double x2) {
    /* Create array of n1 linearly spaced values, plus n2 logarithmically
     * values. */
    assert(x1 >= 0 && x2 >= x1);
    int n = n1 + n2;
    vector<double> X(n);
    for(int i = 0; i < n1; i++)
        X[i] = i * x1 / n1;
    for(int j = 0; j < n2; j++)
        X[n1 + j] = x1 * exp(j*log(x2/x1)/(n2-1));
    return X;
}


QSpline::QSpline() {
    x1 = x2 = y2 = 0;
}

QSpline::QSpline(int n1, int n2, const double x[], const double y[]) {
    int n = n1 + n2;
    assert(n >= 2 && n1 >= 0 && n2 >= 0);
    x1 = x[n1];
    x2 = x[n-1];
    y2 = y[n-1];

    /* Linear section */
    F1 = LinearSpline(n1+1, &x[0], &y[0]);

    /* Logarithmic section: use u = log(x) for indexing */
    vector<double> u(n2);
    for(int j = 0; j < n2; j++)
        u[j] = log(x[n1 + j]);
    F2 = LinearSpline(n2, &u[0], &y[n1]);
}

QSpline::QSpline(int n1, int n2, const vector<double>& x, const vector<double>& y) {
    int n = (int) x.size();
    assert(n >= 2 && n1 >= 0 && n2 >= 0 && n1 + n2 == n && (int) y.size() >= n);
    x1 = x[n1];
    x2 = x[n-1];
    y2 = y[n-1];

    /* Linear section */
    F1 = LinearSpline(n1+1, &x[0], &y[0]);

    /* Logarithmic section: use u = log(x) for indexing */
    vector<double> u(n2);
    for(int j = 0; j < n2; j++)
        u[j] = log(x[n1 + j]);
    F2 = LinearSpline(n2, &u[0], &y[n1]);
}

QSpline::~QSpline() {
}

double QSpline::operator()(double x) const {
    if(x < 0)
        return 0;
    else if(x < x1)
        return F1(x);
    else if(x < x2)
        return F2(log(x));
    else
        return y2;
}
