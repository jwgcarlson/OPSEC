#include <cassert>
#include <cmath>
#include <cstring>

#include "TestModel.h"

#define DEFINE_TEST_XI(Name, formula) \
struct Name : public XiFuncImpl { \
    double A, b; \
    Name(double A_, double b_) : A(A_), b(b_) {} \
    double xi(double r, double, double) const { \
        return formula; \
    } \
};

DEFINE_TEST_XI(ExponentialXi, A*exp(-b*r))
DEFINE_TEST_XI(GaussianXi, A*exp(-b*r*r))
DEFINE_TEST_XI(DGaussianXi, A*r*exp(-b*r*r))


TestModel::TestModel(Config cfg) {
    if(cfg_has_key(cfg, "form")) {
        const char* s = cfg_get(cfg, "form");
        if(strcmp(s, "exponential") == 0)
            form = 0;
        else if(strcmp(s, "gaussian") == 0)
            form = 1;
        else {
            fprintf(stderr, "TestModel: unrecognized functional form: '%s'\n", s);
            fprintf(stderr, "           (defaulting to 'exponential')\n");
            form = 0;
        }
    }
    A = cfg_has_key(cfg, "A") ? cfg_get_double(cfg, "A") : 1.0;
    b = cfg_has_key(cfg, "b") ? cfg_get_double(cfg, "b") : 0.01;

#if 0
    const char* s = cfg_has_key(cfg, "form") ?  cfg_get(cfg, "form") : "exponential";
    if(strcmp(s, "exponential") == 0) {
        double A = cfg_has_key(cfg, "A") ? cfg_get_double(cfg, "A") : 1.0;
        double b = cfg_has_key(cfg, "b") ? cfg_get_double(cfg, "b") : 0.01;
        xi = XiFunc(new ExponentialXi(A, b));
    }
    else if(strcmp(s, "gaussian") == 0) {
        double A = cfg_has_key(cfg, "A") ? cfg_get_double(cfg, "A") : 1.0;
        double b = cfg_has_key(cfg, "b") ? cfg_get_double(cfg, "b") : 0.01;
        xi = XiFunc(new GaussianXi(A, b));
    }
    else {
        fprintf(stderr, "TestModel: Unrecognized functional form: '%s'\n", s);
        fprintf(stderr, "           Defaulting to 'exponential'.\n");
        xi = XiFunc(new ExponentialXi(1.0, 0.01));
    }
#endif
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
