#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "RealSplineModel.h"
#include "opsec.h"

using namespace std;

/* Implementation of a real-space (homogeneous and isotropic) correlation
 * function, splined as xi(r) vs. log(r). */
struct RealSpaceXi : public XiFuncImpl {
    Spline xir;
    double r0, xi0;

    RealSpaceXi(Spline xir_, double r0_, double xi0_) {
        xir = xir_;
        r0 = r0_;
        xi0 = xi0_;
    }

    double xi(double r, double, double) const {
        return (r <= r0) ? xi0 : xir(log(r));
    }
};


RealSplineModel::RealSplineModel(Config cfg) {
    /* Read prior power spectrum from file */
    if(!cfg_has_key(cfg, "pkfile")) {
        fprintf(stderr, "RealSplineModel: must set pkfile in config file\n");
        return;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    pk = LinearSpline(pkfile);

    /* Determine knot locations */
    if(cfg_has_key(cfg, "knots")) {
        /* Read knots directly */
        char* knotstr = strdup(cfg_get(cfg, "knots"));
        char* saveptr;
        char* s = strtok_r(knotstr, " \t", &saveptr);
        while(s != NULL) {
            double k = atof(s);
            knots.push_back(k);
            s = strtok_r(NULL, " \t", &saveptr);
        }
        free(knotstr);
    }
    else if(cfg_has_keys(cfg, "Nknots,knotmin,knotmax", ",")) {
        /* Use regularly spaced knots */
        int Nknots = cfg_get_int(cfg, "Nknots");
        double kmin = cfg_get_double(cfg, "knotmin");
        double kmax = cfg_get_double(cfg, "knotmax");
        /* (Assume linearly spaced knots for now) */
        for(int i = 0; i < Nknots; i++)
            knots.push_back(kmin + i*(kmax - kmin)/(Nknots-1));
    }
    else {
        fprintf(stderr, "RealSplineModel: k knots not specified in config file\n");
        return;
    }

    Nknots = (int) knots.size();

    /* Check knot locations for sanity (TODO) */
}

RealSplineModel::~RealSplineModel() {
}

int RealSplineModel::NumParams() {
    return Nknots;
}

/* Spherical Bessel function $j_0(x)$ */
inline static double sj0(double x) {
    if(fabs(x) < 1e-4) {
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x2*x4;
        return 1 - x2/6 + x4/120 - x6/5040;
    }
    else
        return sin(x)/x;
}

XiFunc RealSplineModel::GetXi() {
    const int Nr = 1024;        // number of points to evaluate xi(r)
    const double rmin = 1e-2;
    const double rmax = 1e4;
    const int Nk = 8192;        // number of _subintervals_, not number of k points; must be even
    const double kmin = 0.;
    const double kmax = 10.;
    double k, dk;
    vector<double> r(Nr), logr(Nr), xi(Nr);
    int mult;

    dk = (kmax - kmin)/Nk;

    for(int i = 0; i < Nr; i++) {
        r[i] = rmin * pow(rmax/rmin, i/(Nr-1.));
        logr[i] = log(r[i]);
        xi[i] = 0;
    }

    /* Integrate $P(k) k^2 j_0(kr) dk$ from kmin to kmax using Simpson's rule */
    for(int j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        for(int i = 0; i < Nr; i++)
            xi[i] += mult * k*k*sj0(k*r[i]) * pk(k);
    }
    for(int i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);

    return XiFunc(new RealSpaceXi(LinearSpline(Nr, &logr[0], &xi[0]), r[0], xi[0]));
}

static inline double tent(double k, double kmin, double kmid, double kmax) {
    int w = (k <= kmid);
    return w * (k-kmin)/(kmid-kmin) + (1-w) * (kmax-k)/(kmax-kmid);
}

XiFunc RealSplineModel::GetXiDeriv(int n) {
    const int Nr = 1024;        // number of points to evaluate xi,m(r)
    const double rmin = 1e-2;
    const double rmax = 1e4;
    const int Nk = 8192;        // number of _subintervals_, not number of k points; must be even
    double k, kmin, kmid, kmax, dk;
    vector<double> r(Nr), logr(Nr), xi(Nr);
    int mult;

    kmin = (n == 0) ? 0. : knots[n-1];
    kmid = knots[n];
    kmax = (n == Nknots-1) ? kmid + (kmid-knots[n-1]) : knots[n+1];
    dk = (kmax - kmin)/Nk;

    for(int i = 0; i < Nr; i++) {
        r[i] = rmin * pow(rmax/rmin, i/(Nr-1.));
        logr[i] = log(r[i]);
        xi[i] = 0;
    }

    /* Integrate $T_n(k) k^2 j_0(kr) dk$ over the interval [kmin,kmax] using
     * Simpson's rule, where $T_n(k)$ is a tent function */
    for(int j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        for(int i = 0; i < Nr; i++)
            xi[i] += mult * k*k*sj0(k*r[i]) * tent(k, kmin, kmid, kmax);
    }
    for(int i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);

    return XiFunc(new RealSpaceXi(LinearSpline(Nr, &logr[0], &xi[0]), r[0], xi[0]));
}
