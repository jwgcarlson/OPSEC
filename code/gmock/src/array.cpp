#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "array.h"

using std::vector;


array::array() {
    n[0] = n[1] = n[2] = 0;
}

array::array(size_type n0) : vector<double>(n0) {
    n[0] = n0;
    n[1] = n[2] = 1;
}

array::array(size_type n0, size_type n1) : vector<double>(n0*n1) {
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

array::array(size_type n0, size_type n1, size_type n2) : vector<double>(n0*n1*n2) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

array::array(size_type n0, const double* v) : vector<double>(v, v+n0) {
    n[0] = n0;
    n[1] = n[2] = 1;
}

array::array(size_type n0, size_type n1, const double* v) : vector<double>(v, v+n0*n1) {
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

array::array(size_type n0, size_type n1, size_type n2, const double* v) : vector<double>(v, v+n0*n1*n2) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

array::array(const array& v) : vector<double>(v) {
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
}

array& array::operator=(const array& v) {
    vector<double>::operator=(v);
    n[0] = v.n[0];
    n[1] = v.n[1];
    n[2] = v.n[2];
    return *this;
}

array::~array() {
}

array::array(size_type n0, size_type n1, size_type n2, double r) : vector<double>(n0*n1*n2, r) {
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

array array::zeros(size_type n0) {
    return array(n0, 1, 1, 0.);
}

array array::zeros(size_type n0, size_type n1) {
    return array(n0, n1, 1, 0.);
}

array array::zeros(size_type n0, size_type n1, size_type n2) {
    return array(n0, n1, n2, 0.);
}

array array::ones(size_type n0) {
    return array(n0, 1, 1, 1.);
}

array array::ones(size_type n0, size_type n1) {
    return array(n0, n1, 1, 1.);
}

array array::ones(size_type n0, size_type n1, size_type n2) {
    return array(n0, n1, n2, 1.);
}

array::array(size_type n0, double min, double max, Spacing mode) : vector<double>(n0) {
    n[0] = n0;
    n[1] = n[2] = 1;

    assert(n0 >= 2);
    if(mode == LinearSpacing) {
        for(size_type i = 0; i < n0; i++)
            (*this)[i] = min + i*(max - min)/(n0-1);
    }
    else if(mode == LogarithmicSpacing) {
        assert(min > 0 && max > 0);
        double logmin = log(min);
        double logmax = log(max);
        for(size_type i = 0; i < n0; i++)
            (*this)[i] = exp(logmin + i*(logmax - logmin)/(n0-1));
    }

    /* Guarantee min and max elements are exact */
    (*this)[0] = min;
    (*this)[n0-1] = max;
}

array array::linspace(double min, double max, size_type n0) {
    return array(n0, min, max, LinearSpacing);
}

array array::logspace(double min, double max, size_type n0) {
    return array(n0, min, max, LogarithmicSpacing);
}

void array::resize(size_type n0) {
    vector<double>::resize(n0);
    n[0] = n0;
    n[1] = n[2] = 1;
}

void array::resize(size_type n0, size_type n1) {
    vector<double>::resize(n0*n1);
    n[0] = n0;
    n[1] = n1;
    n[2] = 1;
}

void array::resize(size_type n0, size_type n1, size_type n2) {
    vector<double>::resize(n0*n1*n2);
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
}

//void array::transpose(int axis0, int axis1) {
//}

void array::getshape(size_type* n0, size_type* n1, size_type* n2) const {
    if(n0) *n0 = n[0];
    if(n1) *n1 = n[1];
    if(n2) *n2 = n[2];
}

double array::min() const {
    size_type N = size();
    if(N == 0)
        return 0.;
    double r = 1e100;
    for(size_type i = 0; i < N; i++)
        r = fmin(r, (*this)[i]);
    return r;
}

double array::max() const {
    size_type N = size();
    if(N == 0)
        return 0.;
    double r = -1e100;
    for(size_type i = 0; i < N; i++)
        r = fmax(r, (*this)[i]);
    return r;
}

array array::abs() const {
    size_type N = size();
    if(N == 0)
        return array();
    array r(N);
    for(size_type i = 0; i < N; i++)
        r[i] = fabs((*this)[i]);
    return r;
}

#if 0
void array::setdata(size_type n, double r) {
    resize(n, r);
}

void array::setdata(size_type n, const double* v) {
    resize(n);
    for(iterator it = begin(); it != end(); it++)
        *it = *v++;
}
#endif

array& array::operator+=(const array& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] += v[i];
    return *this;
}

array& array::operator-=(const array& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] -= v[i];
    return *this;
}

array& array::operator*=(const array& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] *= v[i];
    return *this;
}

array& array::operator/=(const array& v) {
    size_type N = size();
    assert(v.size() == N);
    for(size_type i = 0; i < N; i++)
        (*this)[i] /= v[i];
    return *this;
}

array& array::operator+=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] += s;
    return *this;
}

array& array::operator-=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] -= s;
    return *this;
}

array& array::operator*=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] *= s;
    return *this;
}

array& array::operator/=(double s) {
    size_type N = size();
    for(size_type i = 0; i < N; i++)
        (*this)[i] /= s;
    return *this;
}

array operator-(const array& v) {
    array::size_type N = v.size();
    array w(v);
    for(array::size_type i = 0; i < N; i++)
        w[i] = -w[i];
    return w;
}

array operator+(const array& u, const array& v) {
    assert(u.size() == v.size());
    array w = u;
    w += v;
    return w;
}

array operator-(const array& u, const array& v) {
    assert(u.size() == v.size());
    array w = u;
    w -= v;
    return w;
}

array operator*(const array& u, const array& v) {
    assert(u.size() == v.size());
    array w = u;
    w *= v;
    return w;
}

array operator/(const array& u, const array& v) {
    assert(u.size() == v.size());
    array w = u;
    w /= v;
    return w;
}

array operator+(const array& u, double s) {
    array w = u;
    w += s;
    return w;
}

array operator+(double s, const array& u) {
    array w = u;
    w += s;
    return w;
}

array operator-(const array& u, double s) {
    array w = u;
    w -= s;
    return w;
}

array operator*(const array& u, double s) {
    array w = u;
    w *= s;
    return w;
}

array operator*(double s, const array& u) {
    array w = u;
    w *= s;
    return w;
}

array operator/(const array& u, double s) {
    array w = u;
    w /= s;
    return w;
}

double min(const array& v) {
    return v.min();
}

double max(const array& v) {
    return v.max();
}

array abs(const array& v) {
    return v.abs();
}
