#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "XiModel.h"
#include "opsec.h"

using namespace std;

/* Implementation of a real-space (homogeneous and isotropic) correlation
 * function */
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


XiModel::XiModel(Config cfg) {

    /* Read prior power spectrum from file */
    if(!cfg_has_key(cfg, "XiModel.pkfile")) {
        fprintf(stderr, "XiModel: must set XiModel.pkfile in config file\n");
        return;
    }
    const char* pkfile = cfg_get(cfg, "XiModel.pkfile");
    pk = CubicSpline(pkfile);

    /* Determine r bands */
    if(cfg_has_key(cfg, "XiModel.rbands")) {
        /* Read bands directly */
        char* rbandstr = strdup(cfg_get(cfg, "XiModel.rbands"));
        char* saveptr;
        char* s = strtok_r(rbandstr, " \t", &saveptr);
        while(s != NULL) {
            double r = atof(s);
            rbands.push_back(r);
            s = strtok_r(NULL, " \t", &saveptr);
        }
        free(rbandstr);
    }
    else if(cfg_has_keys(cfg, "XiModel.Nbands,XiModel.rbandmin,XiModel.rbandmax", ",")) {
        /* Use regularly spaced bands */
        int Nbands = cfg_get_int(cfg, "XiModel.Nbands");
        double rmin = cfg_get_double(cfg, "XiModel.rbandmin");
        double rmax = cfg_get_double(cfg, "XiModel.rbandmax");
        /* (Assume linearly spaced bands for now) */
        for(int i = 0; i <= Nbands; i++)
            rbands.push_back(rmin + i*(rmax - rmin)/Nbands);
    }
    else {
        fprintf(stderr, "XiModel: r bands not specified in config file\n");
        return;
    }

    /* Verify the r bands are sane (TODO) */
}

XiModel::~XiModel() {
}

int XiModel::NumParams() {
    return (int) rbands.size() - 1;
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

XiFunc XiModel::GetXi() {
    const int Nr = 1024;        // number of points to evaluate xi(r)
    const double rmin = 1e-2;
    const double rmax = 1e4;
    const int Nk = 8192;     // number of _subintervals_, not number of k points; must be even
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

    /* Integrate $k^2 j_0(kr) dk$ over the interval $[k_n,k_{n+1}]$ using Simpson's rule */
    for(int j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        for(int i = 0; i < Nr; i++)
            xi[i] += mult * pk(k) * k*k*sj0(k*r[i]);
    }
    for(int i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);

    return XiFunc(new RealSpaceXi(LinearSpline(Nr, &logr[0], &xi[0]), r[0], xi[0]));
}

XiFunc XiModel::GetXiDeriv(int n) {
    const int Nr = 16384;        // number of points to evaluate \xi,m(r)
    const double rmin = 1e-2;
    const double rmax = 1e4;

    vector<double> r(Nr), logr(Nr), xi(Nr);
    for(int i = 0; i < Nr; i++) {
        r[i] = rmin * pow(rmax/rmin, i/(Nr-1.));
        logr[i] = log(r[i]);
        xi[i] = (rbands[n] <= r[i] && r[i] < rbands[n+1]);    // either 1 or 0
    }

    return XiFunc(new RealSpaceXi(LinearSpline(Nr, &logr[0], &xi[0]), r[0], xi[0]));
}
