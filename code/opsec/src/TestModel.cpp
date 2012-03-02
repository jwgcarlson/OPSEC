#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstring>

#include "SeparationFunc.h"
#include "TestModel.h"
#include "XiFunc.h"

#define DEFINE_TEST_XI(Name, formula) \
struct Name : public XiFuncImpl { \
    double A, b; \
    Name(double A_, double b_) : A(A_), b(b_) {} \
    double xi(const Point& p1, const Point& p2, const SeparationFunc& sep) { \
        double r = sep.r(p1, p2); \
        return formula; \
    } \
};

DEFINE_TEST_XI(ExponentialXi, A*exp(-b*r))
DEFINE_TEST_XI(GaussianXi, A*exp(-b*r*r))
DEFINE_TEST_XI(DGaussianXi, A*r*exp(-b*r*r))


TestModel::TestModel(int form_, double A_, double b_) {
    form = form_;
    A = A_;
    b = b_;
}

TestModel::TestModel(Config cfg) {
    form = cfg_get_enum(cfg, "form", "exponential", 0,
                                     "gaussian", 1,
                                     "", -1);
    if(form == -1) {
        fprintf(stderr, "TestModel: unrecognized functional form: '%s'\n", cfg_get(cfg, "form"));
        fprintf(stderr, "           (defaulting to 'exponential')\n");
        form = 0;
    }
    A = cfg_has_key(cfg, "A") ? cfg_get_double(cfg, "A") : 1.0;
    b = cfg_has_key(cfg, "b") ? cfg_get_double(cfg, "b") : 0.01;
}

TestModel::~TestModel() {
}

int TestModel::NumParams() {
    return 2;
}

double TestModel::GetParam(int n) {
    assert(0 <= n && n < 2);
    switch(n) {
    case 0:
        return A;
    case 1:
        return b;
    }
}

XiFunc TestModel::GetXi() {
    if(form == 0)
        return XiFunc(new ExponentialXi(A, b));
    else
        return XiFunc(new GaussianXi(A, b));
}

XiFunc TestModel::GetXiDeriv(int n) {
    assert(0 <= n && n < 2);
    switch(n) {
    case 0:
        if(form == 0)
            return XiFunc(new ExponentialXi(1., b));
        else
            return XiFunc(new GaussianXi(1., b));
    case 1:
        if(form == 0)
            return XiFunc(new ExponentialXi(-A*b, b));
        else
            return XiFunc(new DGaussianXi(-A*b, b));
    }
}
