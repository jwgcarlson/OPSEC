#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <vector>
using std::vector;

#include "KaiserModel.h"
#include "opsec.h"

/* Implementation of a redshift-space correlation function in the
 * plane-parallel approximation. */
struct PlaneParallelXi : public XiFuncImpl {
    Spline xi0, xi2, xi4;
    double smin;

    PlaneParallelXi(Spline xi0_, Spline xi2_, Spline xi4_, double smin_) {
        xi0 = xi0_;
        xi2 = xi2_;
        xi4 = xi4_;
        smin = smin_;
    }

    /* See src/Model.h or doc/geometry.tex for details. */
    double xi(double s, double a, double b) const {
        /* We avoid inverse trigonometric functions, as they tend to be slow
         * and this method needs to be as fast as possible. */
        double cos2theta = fmin(1, (a*a + b*b - s*s)/(2*a*b));
        double costheta = sqrt(0.5*(1 + cos2theta));
        double sintheta = sqrt(0.5*(1 - cos2theta));
        double cosA = fmin(1, (s*s + b*b - a*a)/(2*s*b));
        double sinA = sqrt(1 - cosA*cosA);
        double mu = cosA*costheta - sinA*sintheta, mu2 = mu*mu, mu4 = mu2*mu2;
        double L0 = 1., L2 = (3*mu2 - 1)/2., L4 = (35*mu4 - 30*mu2 + 3)/8.;

        double logs = log(fmax(s,smin));
        return L0*xi0(logs) + L2*xi2(logs) + L4*xi4(logs);
    }
};


KaiserModel::KaiserModel(Config cfg) {
    /* Read prior power spectrum from file */
    if(!cfg_has_key(cfg, "pkfile")) {
        fprintf(stderr, "KaiserModel: must set config option 'pkfile'\n");
        return;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    pk = CubicSpline(pkfile);

    /* Get bias and growth rate */
    if(!cfg_has_keys(cfg, "f,b", ",")) {
        fprintf(stderr, "KaiserModel: must set config options 'f' and 'b'\n");
        return;
    }
    f = cfg_get_double(cfg, "f");
    b = cfg_get_double(cfg, "b");

    /* Determine bands */
    Band b;
    if(cfg_has_key(cfg, "bands")) {
        /* Read bands directly */
        vector<double> kvals;
        char* bandstr = strdup(cfg_get(cfg, "bands"));
        char* saveptr;
        char* s = strtok_r(bandstr, " \t", &saveptr);
        while(s != NULL) {
            double k = atof(s);
            kvals.push_back(k);
            s = strtok_r(NULL, " \t", &saveptr);
        }
        free(bandstr);

        if((kvals.size() % 2) != 0) {
            fprintf(stderr, "KaiserModel: odd number of kvals\n");
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
        fprintf(stderr, "KaiserModel: k-bands not specified in configuration\n");
        return;
    }
}

KaiserModel::~KaiserModel() {
}

int KaiserModel::NumParams() {
    return 3*Nbands;
}

/* Spherical Bessel functions $j_\ell(x)$ */
inline static double sj0(double x) {
    if(fabs(x) < 1e-4) {
        double x2 = x*x, x4 = x2*x2;
        return 1 - x2/6 + x4/120;
    }
    else
        return sin(x)/x;
}
inline static double sj2(double x) {
    if(fabs(x) < 1e-3) {
        double x2 = x*x, x4 = x2*x2, x6 = x2*x4;
        return x2/15 - x4/210 + x6/7560;
    }
    else {
        double x2 = x*x, x3 = x2*x;
        return ((3 - x2)*sin(x) - 3*x*cos(x))/x3;
    }
}
inline static double sj4(double x) {
    if(fabs(x) < 1e-2) {
        double x2 = x*x, x4 = x2*x2, x6 = x2*x4, x8 = x4*x4;
        return x4/945 - x6/20790 + x8/1081080;
    }
    else {
        double x2 = x*x, x4 = x2*x2, x5 = x4*x;
        return ((x4 - 45*x2 + 105)*sin(x) + x*(10*x2 - 105)*cos(x))/x5;
    }
}

XiFunc KaiserModel::GetXi() {
    const int Ns = 1024;        // number of points to evaluate \xi_\ell(s) before splining
    const double smin = 1e-2;
    const double smax = 1e4;
    const int Nk = 8192;        // number of _subintervals_, not number of k points; must be even
    const double kmin = 0.;
    const double kmax = 1.;

    vector<double> s(Ns), logs(Ns), xi0(Ns, 0.), xi2(Ns, 0.), xi4(Ns, 0.);

    /* Logarithmically spaced s values */
    double logsmin = log(smin), logsmax = log(smax);
    for(int i = 0; i < Ns; i++) {
        logs[i] = logsmin + i*(logsmax - logsmin)/(Ns-1);
        s[i] = exp(logs[i]);
    }

    /* Integrate $k^2 P_\ell(k) j_\ell(ks)$ over $[kmin,kmax]$ using Simpson's rule */
    double k, dk = (kmax-kmin)/Nk, mult;
    for(int j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        mult *= k*k * pk(k);
        for(int i = 0; i < Ns; i++) {
            xi0[i] += mult * sj0(k*s[i]);
            xi2[i] += mult * sj2(k*s[i]);
            xi4[i] += mult * sj4(k*s[i]);
        }
    }

    /* Proportionality factors in $P_\ell(k) = A_\ell P(k)$, in the Kaiser limit */
    double beta = f/b;
    double A0 = 1 + (2/3.)*beta + (1/5.)*beta*beta;
    double A2 = (4/3.)*beta + (4/7.)*beta*beta;
    double A4 = (8/35.)*beta*beta;
    for(int i = 0; i < Ns; i++) {
        xi0[i] *=  A0 * (dk/3) / (2*M_PI*M_PI);
        xi2[i] *= -A2 * (dk/3) / (2*M_PI*M_PI); // minus because of $i^\ell$
        xi4[i] *=  A4 * (dk/3) / (2*M_PI*M_PI);
    }

    return XiFunc(new PlaneParallelXi(LinearSpline(logs, xi0),
                                      LinearSpline(logs, xi2),
                                      LinearSpline(logs, xi4),
                                      smin));
}

XiFunc KaiserModel::GetXiDeriv(int n) {
    assert(0 <= n && n < NumParams());

    const int Ns = 1024;        // number of points to evaluate \xi_\ell(s) before splining
    const double smin = 1e-2;
    const double smax = 1e4;
    const int Nk = 8192;        // number of _subintervals_, not number of k points; must be even

    int ell = 2*(n/Nbands);     // either 0, 2, or 4
    int m = n % Nbands;         // 0 <= m < Nbands
    double kmin = bands[m].min;
    double kmax = bands[m].max;

    vector<double> s(Ns), logs(Ns), xi0(Ns, 0.), xi2(Ns, 0.), xi4(Ns, 0.);

    /* Logarithmically spaced s values */
    double logsmin = log(smin), logsmax = log(smax);
    for(int i = 0; i < Ns; i++) {
        logs[i] = logsmin + i*(logsmax - logsmin)/(Ns-1);
        s[i] = exp(logs[i]);
    }

    /* Integrate $k^2 P_\ell(k) j_\ell(ks)$ over $[kmin,kmax]$ using Simpson's rule */
    double k, dk = (kmax-kmin)/Nk, mult;
    for(int j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        mult *= k*k;
        for(int i = 0; i < Ns; i++) {
            /* Only want \xi_\ell to be nonzero: this is a cheap way to do it */
            xi0[i] += (ell == 0) * mult * sj0(k*s[i]);
            xi2[i] += (ell == 2) * mult * sj2(k*s[i]);
            xi4[i] += (ell == 4) * mult * sj4(k*s[i]);
        }
    }

    /* Proportionality factors in $P_\ell(k) = A_\ell P(k)$, in the Kaiser limit */
    double beta = f/b;
    double A0 = 1 + (2/3.)*beta + (1/5.)*beta*beta;
    double A2 = (4/3.)*beta + (4/7.)*beta*beta;
    double A4 = (8/35.)*beta*beta;
    for(int i = 0; i < Ns; i++) {
        xi0[i] *=  A0 * (dk/3) / (2*M_PI*M_PI);
        xi2[i] *= -A2 * (dk/3) / (2*M_PI*M_PI); // minus because of $i^\ell$
        xi4[i] *=  A4 * (dk/3) / (2*M_PI*M_PI);
    }

    return XiFunc(new PlaneParallelXi(LinearSpline(logs, xi0),
                                      LinearSpline(logs, xi2),
                                      LinearSpline(logs, xi4),
                                      smin));
}
