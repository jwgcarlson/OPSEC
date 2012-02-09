#include <cmath>
#include <cstdio>
#include <cstring>

#include <vector>
using std::vector;

#include "Model.h"

/* Should these somehow use "#pragma omp critical"'s to prevent threading issues? */

XiFunc::XiFunc(XiFuncImpl* impl_) {
    impl = impl_;
    if(impl)
        impl->refcount++;
}

XiFunc::XiFunc(const XiFunc& xi) {
    impl = xi.impl;
    if(impl)
        impl->refcount++;
}

XiFunc& XiFunc::operator=(const XiFunc& xi) {
    if(impl != xi.impl) {
        if(impl && --impl->refcount <= 0)
            delete impl;
        impl = xi.impl;
        if(impl)
            impl->refcount++;
    }
    return *this;
}

XiFunc::~XiFunc() {
    if(impl && --impl->refcount <= 0)
        delete impl;
}

double XiFunc::operator()(double r, double a, double b) const {
#ifdef DEBUG
    return impl ? impl->xi(r, a, b) : 0.0;
#else
    return impl->xi(r, a, b);
#endif
}


Model::Model() {
}

Model::~Model() {
}


//#include "KaiserModel.h"
#include "RealModel.h"
//#include "RealSplineModel.h"
#include "TestModel.h"
//#include "XiModel.h"

Model* InitializeModel(Config cfg) {
    const char* model = cfg_get(cfg, "model");
    Model* m = NULL;
    Config subcfg = NULL;

#define HANDLE_MODEL(NAME) \
    if(strcmp(model, #NAME) == 0) { \
        subcfg = cfg_new_sub(cfg, #NAME ".", 1); \
        m = new NAME(subcfg); \
    }

    HANDLE_MODEL(TestModel)
//    else HANDLE_MODEL(KaiserModel)
    else HANDLE_MODEL(RealModel)
//    else HANDLE_MODEL(RealSplineModel)
//    else HANDLE_MODEL(SzapudiModel)
//    else HANDLE_MODEL(XiModel)
    else {
        fprintf(stderr, "InitializeModel: unrecognized model '%s'\n", model);
        fflush(stderr);
    }
#undef HANDLE_MODEL

    cfg_destroy(subcfg);
    return m;
}


/* Spherical Bessel functions $j_n(x)$ for 0 <= n <= 4. */
static double sj0(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? 1 - x2*(1/6. + x2*(1/120. - x2*(1/5040.)))
                       : sin(x)/x;
}
static double sj1(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x*(1/3. - x2*(1/30. + x2*(1/840.)))
                       : -(x*cos(x) - sin(x))/x2;
}
static double sj2(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x2*(1/15. - x2*(1/210. + x2*(1/7560.)))
                       : -(3.*x*cos(x) + (x2 - 3.)*sin(x))/(x*x2);
}
static double sj3(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x*x2*(1/105. - x2*(1/1890.))
                       : (x*(x2 - 15.)*cos(x) - (6.*x2 - 15.)*sin(x))/(x2*x2);
}
static double sj4(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? x2*x2*(1/945. - x2*(1/20790.))
                       : (x*(10.*x2 - 105.)*cos(x) + (x2*(x2 - 45.) + 105.)*sin(x))/(x*x2*x2);
}

/* Compute the quantity
 *   \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
 * using direct Simpson integration. */
void ComputeXiLM(int l, int m, Spline P, int Nr, const double r[], double xi[],
                 int Nk, double kmin, double kmax, double kdamp)
{
    const double dk = (kmax - kmin)/Nk;

    /* Choose appropriate spherical Bessel function */
    double (*sj)(double x);
    if(l == 0)      sj = sj0;
    else if(l == 1) sj = sj1;
    else if(l == 2) sj = sj2;
    else if(l == 3) sj = sj3;
    else if(l == 4) sj = sj4;
    else {
        fprintf(stderr, "ComputeXiLM: l = %d not supported\n", l);
        return;
    }

    vector<double> k(Nk+1);
    vector<double> mult(Nk+1);

    /* Pre-compute k-dependent factors */
    #pragma omp parallel for
    for(int j = 0; j <= Nk; j++) {
        k[j] = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult[j] = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        /* Other k-dependent factors */
        mult[j] *= P(k[j]) * pow(k[j], m) * exp(-k[j]*k[j]/(kdamp*kdamp));
        /* Normalization factors */
        mult[j] *= (dk/3) / (2*M_PI*M_PI);
    }

    /* Integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
    #pragma omp parallel for
    for(int i = 0; i < Nr; i++) {
        xi[i] = 0;
        for(int j = 0; j <= Nk; j++)
            xi[i] += mult[j] * sj(k[j]*r[i]);
    }
}
