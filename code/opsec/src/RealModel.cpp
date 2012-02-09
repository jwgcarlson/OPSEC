#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "RealModel.h"
#include "QSpline.h"
#include "opsec.h"

/* Implementation of a real-space (homogeneous and isotropic) correlation
 * function */
struct RealSpaceXi : public XiFuncImpl {
    QSpline F;

    RealSpaceXi(QSpline F_) {
        F = F_;
    }

    double xi(double r, double, double) const {
        return F(r);
    }
};


RealModel::RealModel(Config cfg) {
    /* Read prior power spectrum from file */
    if(!cfg_has_key(cfg, "pkfile")) {
        fprintf(stderr, "RealModel: must set config option 'pkfile'\n");
        return;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    pk = CubicSpline(pkfile);

    /* Determine bands */
    Band b;
    if(cfg_has_key(cfg, "bands")) {
        /* Read bands directly */
//        vector<double> kvals;
//        char* bandstr = strdup(cfg_get(cfg, "bands"));
//        char* saveptr;
//        char* s = strtok_r(bandstr, " \t", &saveptr);
//        while(s != NULL) {
//            double k = atof(s);
//            kvals.push_back(k);
//            s = strtok_r(NULL, " \t", &saveptr);
//        }
//        free(bandstr);
        vector<double> kvals(100);
        cfg_get_array_double(cfg, "bands", 100, &kvals[0]);

        if((kvals.size() % 2) != 0) {
            fprintf(stderr, "RealModel: odd number of kvals\n");
            kvals.pop_back();
        }
        for(int n = 0; n < (int)kvals.size()/2; n++) {
            b.min = kvals[2*n];
            b.max = kvals[2*n+1];
            bands.push_back(b);
        }
        Nbands = bands.size();
    }
    else if(cfg_has_keys(cfg, "Nbands,kmin,kmax", ",")) {
        /* Use regularly spaced bands */
        Nbands = cfg_get_int(cfg, "Nbands");
        double kmin = cfg_get_double(cfg, "kmin");
        double kmax = cfg_get_double(cfg, "kmax");
        /* (Assume linearly spaced bands for now) */
        for(int n = 0; n < Nbands; n++) {
            b.min = kmin + n*(kmax - kmin)/Nbands;
            b.max = kmin + (n+1)*(kmax - kmin)/Nbands;
            bands.push_back(b);
        }
    }
    else {
        fprintf(stderr, "RealModel: k-bands not specified in configuration\n");
        return;
    }
}

RealModel::~RealModel() {
}

int RealModel::NumParams() {
    return Nbands;
}

double RealModel::GetParam(int n) {
    assert(0 <= n && n < Nbands);
    const int Nk = 256;
    double k, kmin = bands[n].min, kmax = bands[n].max;
    /* Compute average value of P(k) over the interval [kmin,kmax], i.e. a "band power" */
    double sum = 0;
    for(int i = 0; i < Nk; i++) {
        k = kmin + (i + 0.5)*(kmax - kmin)/Nk;
        sum += pk(k);
    }
    return sum/Nk;
}

XiFunc RealModel::GetXi() {
    const int n1 = 512;
    const int n2 = 256;
    const double r1 = 200;
    const double r2 = 20000;

    vector<double> r = QSpline::MakeArray(n1, n2, r1, r2);
    vector<double> xi(n1+n2);
    ComputeXiLM(0, 2, pk, n1+n2, &r[0], &xi[0]);
    return XiFunc(new RealSpaceXi(QSpline(n1, n2, r, xi)));
}

XiFunc RealModel::GetXiDeriv(int n) {
    assert(0 <= n && n < Nbands);

    const int Nk = 8192;
    double kmin = bands[n].min;
    double kmax = bands[n].max;

    const int n1 = 512;
    const int n2 = 256;
    const double r1 = 200;
    const double r2 = 20000;

    vector<double> r = QSpline::MakeArray(n1, n2, r1, r2);
    vector<double> xi(n1+n2);
    Spline one = ConstantSpline(1, kmin, kmax);
    ComputeXiLM(0, 2, one, n1+n2, &r[0], &xi[0], Nk, kmin, kmax);
    return XiFunc(new RealSpaceXi(QSpline(n1, n2, r, xi)));
}

#if 0
/* Spherical Bessel function $j_0(x)$ */
inline static double sj0(double x) {
    double x2 = x*x;
    return (x2 < 1e-4) ? 1 - x2*(1/6. + x2*(1/120. - x2*(1/5040.))) : sin(x)/x;
}

XiFunc RealModel::GetXi() {
    const int Nr = 1024;        // number of points to evaluate \xi(r) before splining
    const double rmin = 1e-2;
    const double rmax = 1e4;
    const int Nk = 65536;       // number of _subintervals_, not number of k points; must be even
    const double kmin = 0.;
    const double kmax = 10.;    // this needs to be surprisingly large, to eliminate ringing in the j_0(kr) integral
    int i, j;

    vector<double> r(Nr), logr(Nr), xi(Nr, 0.);

    /* Logarithmically spaced r values */
    double logrmin = log(rmin), logrmax = log(rmax);
    for(i = 0; i < Nr; i++) {
        logr[i] = logrmin + i*(logrmax - logrmin)/(Nr-1);
        r[i] = exp(logr[i]);
    }

    /* Integrate $k^2 P(k) j_0(kr)$ over $[kmin,kmax]$ using Simpson's rule */
    double k, dk = (kmax-kmin)/Nk, mult;
    #pragma omp parallel for private(i,j,k,mult)
    for(j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        mult *= k*k * pk(k);
        for(i = 0; i < Nr; i++)
            #pragma omp atomic
            xi[i] += mult * sj0(k*r[i]);
    }

    #pragma omp parallel for
    for(i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);

    return XiFunc(new RealSpaceXi(LinearSpline(logr, xi), rmin));
}

XiFunc RealModel::GetXiDeriv(int n) {
    assert(0 <= n && n < Nbands);

    const int Nr = 1024;        // number of points to evaluate \xi(r) before splining
    const double rmin = 1e-2;
    const double rmax = 1e4;
    const int Nk = 65536;       // number of _subintervals_, not number of k points; must be even
    double kmin = bands[n].min;
    double kmax = bands[n].max;
    int i, j;

    vector<double> r(Nr), logr(Nr), xi(Nr, 0.);
//    double* xiptr = &xi[0];     // guarantee omp atomic is safe

    /* Logarithmically spaced r values */
    double logrmin = log(rmin), logrmax = log(rmax);
    for(i = 0; i < Nr; i++) {
        logr[i] = logrmin + i*(logrmax - logrmin)/(Nr-1);
        r[i] = exp(logr[i]);
    }

    /* Integrate $k^2 j_0(kr)$ over $[kmin,kmax]$ using Simpson's rule */
    double k, dk = (kmax-kmin)/Nk, mult;
    #pragma omp parallel for private(i,j,k,mult)
    for(j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        mult *= k*k;
        for(i = 0; i < Nr; i++)
            #pragma omp atomic
            xi[i] += mult * sj0(k*r[i]);
//            xiptr[i] += mult * sj0(k*r[i]);
    }

    #pragma omp parallel for
    for(i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);
//        xiptr[i] *= (dk/3) / (2*M_PI*M_PI);

    return XiFunc(new RealSpaceXi(LinearSpline(logr, xi), rmin));
}
#endif
