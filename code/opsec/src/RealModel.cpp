#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "RealModel.h"
#include "QSpline.h"
#include "SeparationFunc.h"
#include "XiFunc.h"
#include "opsec.h"

/* Implementation of a real-space (homogeneous and isotropic) correlation
 * function */
struct RealSpaceXi : public XiFuncImpl {
    QSpline F;

    RealSpaceXi(QSpline F_) {
        F = F_;
    }

    double xi(const Point& p1, const Point& p2, const SeparationFunc& sep) {
        double r = sep.r(p1, p2);
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
