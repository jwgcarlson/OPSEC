#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
using std::vector;

#include "SzapudiModel.h"
#include "opsec.h"

static inline double clamp(double x, double min, double max) {
    if(x < min) return min;
    else if(x > max) return max;
    else return x;
}

/* Compute the Papai-Szapudi linear theory wide-angle redshift-space
 * correlation function. */
static double ComputeSzapudiXi(double r, double r1, double r2, double f,
                               double xi00, double xi02, double xi11, double xi20,
                               double xi22, double xi31, double xi42)
{
    double f2 = f*f;
    double g1 = r1/r;
    double g2 = r2/r;

    double c00 = (1 + 2*f/3. + 2*f2/15.) * xi02;
    double c02 = (-f/2. - 3*f2/14.)*xi22 + f2/28.*xi42;
    double c20 = c02;
    double c22 = f2/15.*xi02 - f2/21.*xi22 + 19*f2/140.*xi42;
    double c01 = ((10*f + 4*f2)*xi11 - f2*xi31)/(5*g1*r);
    double c10 = ((10*f + 4*f2)*xi11 - f2*xi31)/(5*g2*r);
    double c11 = -(4*f2*xi00 - 8*f2*xi20)/(3*g1*g2*r*r);
    double c21 = (2*f2*xi11 - 3*f2*xi31)/(5*g1*r);
    double c12 = (2*f2*xi11 - 3*f2*xi31)/(5*g2*r);

    double d22 = -f2/15.*xi02 + f2/21.*xi20 + 4*f2/35.*xi42;
    double d11 = 4*f2*(xi00 + xi20)/(3*g1*g2*r*r);
    double d21 = -2*f2*(xi11 + xi31)/(5*g1*r);
    double d12 = -2*f2*(xi11 + xi31)/(5*g2*r);

    double cosa = clamp((r*r + r2*r2 - r1*r1)/(2*r*r2), -1, +1);        // \cos(\alpha)
    double cosb = clamp((r*r + r1*r1 - r2*r2)/(2*r*r1), -1, +1);        // \cos(\beta)
    double sina = sqrt(1 - cosa*cosa);
    double sinb = sqrt(1 - cosb*cosb);
    double cos2a = 2*cosa*cosa - 1;
    double cos2b = 2*cosb*cosb - 1;
    double sin2a = 2*sina*cosa;
    double sin2b = 2*sinb*cosb;

    return c00 + c02*cos2b + c20*cos2a + c22*cos2a*cos2b + c01*cosb + c10*cosa
        + c11*cosa*cosb + c21*cos2a*cosb + c12*cosa*cos2b
        + d22*sin2a*sin2b + d11*sina*sinb + d21*sin2a*sinb + d12*sina*sin2b;
}

/* Compute the derivative of the Papai-Szapudi correlation function with
 * respect to the growth rate $f$. */
static double ComputeSzapudiXi_f(double r, double r1, double r2, double f,
                                 double xi00, double xi02, double xi11, double xi20,
                                 double xi22, double xi31, double xi42)
{
    double g1 = r1/r;
    double g2 = r2/r;

    double c00 = (2/3. + 4*f/15.) * xi02;
    double c02 = (-1/2. - 3*f/7.)*xi22 + f/14.*xi42;
    double c20 = c02;
    double c22 = 2*f/15.*xi02 - 2*f/21.*xi22 + 19*f/70.*xi42;
    double c01 = ((10 + 8*f)*xi11 - 2*f*xi31)/(5*g1*r);
    double c10 = ((10 + 8*f)*xi11 - 2*f*xi31)/(5*g2*r);
    double c11 = -(8*f*xi00 - 16*f*xi20)/(3*g1*g2*r*r);
    double c21 = (4*f*xi11 - 6*f*xi31)/(5*g1*r);
    double c12 = (4*f*xi11 - 6*f*xi31)/(5*g2*r);

    double d22 = -2*f/15.*xi02 + 2*f/21.*xi20 + 8*f/35.*xi42;
    double d11 = 8*f*(xi00 + xi20)/(3*g1*g2*r*r);
    double d21 = -4*f*(xi11 + xi31)/(5*g1*r);
    double d12 = -4*f*(xi11 + xi31)/(5*g2*r);

    double cosa = (r*r + r2*r2 - r1*r1)/(2*r*r2);       // \cos(\alpha)
    double cosb = (r*r + r1*r1 - r2*r2)/(2*r*r1);       // \cos(\beta)
    double sina = sqrt(1 - cosa*cosa);
    double sinb = sqrt(1 - sinb*sinb);
    double cos2a = 2*cosa*cosa - 1;
    double cos2b = 2*cosb*cosb - 1;
    double sin2a = 2*sina*cosa;
    double sin2b = 2*sinb*cosb;

    return c00 + c02*cos2b + c20*cos2a + c22*cos2a*cos2b + c01*cosb + c10*cosa
        + c11*cosa*cosb + c21*cos2a*cosb + c12*cosa*cos2b
        + d22*sin2a*sin2b + d11*sina*sinb + d21*sin2a*sinb + d12*sina*sin2b;
}


/* Papai-Szapudi wide-angle redshift space correlation function $\xi(r,a,b)$. */
struct SzapudiXi : public XiFuncImpl {
    /* Interpolated functions $\xi_l^m(r)$, in terms of $x = \log(r)$ */
    Spline xi00, xi02, xi11, xi20, xi22, xi31, xi42;

    /* Logarithmic growth rate $f = d\ln D/d\ln a$ */
    double f;

    /* Return $\partial\xi/\partial f$ instead of $\xi$ */
    bool partial_f;

    SzapudiXi() {
        /* Constructor is empty, structure must be initialized manually */
    }

    double xi(double r, double r1, double r2) const {
        double x = log(r);
        if(partial_f)
            return ComputeSzapudiXi_f(r, r1, r2, f, xi00(x), xi02(x), xi11(x), xi20(x), xi22(x), xi31(x), xi42(x));
        else
            return ComputeSzapudiXi(r, r1, r2, f, xi00(x), xi02(x), xi11(x), xi20(x), xi22(x), xi31(x), xi42(x));
    }
};


SzapudiModel::SzapudiModel(Config cfg) {
    Nbands = 0;

    /* Read prior power spectrum from file */
    if(!cfg_has_key(cfg, "pkfile")) {
        fprintf(stderr, "SzapudiModel: must set pkfile in config file\n");
        return;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    pk = LinearSpline(pkfile);

    /* Get logarithmic growth rate */
    if(!cfg_has_key(cfg, "f")) {
        fprintf(stderr, "SzapudiModel: must set f in config file\n");
        return;
    }
    f = cfg_get_double(cfg, "f");

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

        if((kvals.size() % 2) != 0)
            fprintf(stderr, "SzapudiModel: odd number of kvals\n");
        for(int n = 0; n < kvals.size()/2; n++) {
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
        fprintf(stderr, "SzapudiModel: k-bands not specified in config file\n");
        return;
    }
}

SzapudiModel::~SzapudiModel() {
}

int SzapudiModel::NumParams() {
    /* P(k) bands plus the anisotropy parameter f */
    return Nbands + 1;
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
static void ComputeXiLM(int l, int m, double kmin, double kmax, Spline P, int Nr, const double r[], double xi[]) {
    const int Nk = 65536;     // number of _subintervals_, not number of points; must be even
    const double dk = (kmax - kmin)/Nk;
    int i, j;
    double k, mult;

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

    for(i = 0; i < Nr; i++)
        xi[i] = 0;

    /* Integrate $P(k) k^m j_l(kr) dk$ over the interval $[kmin,kmax]$ using Simpson's rule */
    #pragma omp parallel for private(i,k,mult)
    for(j = 0; j <= Nk; j++) {
        k = kmin + j*dk;
        /* Multiplicative factor for Simpson's rule: either 1, 2, or 4 */
        mult = 2 + 2*(j % 2) - (j == 0) - (j == Nk);
        mult *= P(k) * pow(k, m);
        for(i = 0; i < Nr; i++)
            #pragma omp atomic
            xi[i] += mult * sj(k*r[i]);
    }
    for(i = 0; i < Nr; i++)
        xi[i] *= (dk/3) / (2*M_PI*M_PI);
}

static void FillSzapudiXi(SzapudiXi* szapudi, int Nr, double rmin, double rmax, double kmin, double kmax, Spline P) {
    double logrmin = log(rmin), logrmax = log(rmax);
    vector<double> r(Nr), logr(Nr), xitmp(Nr);
    for(int i = 0; i < Nr; i++) {
        logr[i] = logrmin + i*(logrmax - logrmin)/(Nr-1);
        r[i] = exp(logr[i]);
    }

    /* Pre-compute functions $\xi_l^m(r)$ */
    ComputeXiLM(0, 0, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi00 = LinearSpline(logr, xitmp);
    ComputeXiLM(0, 2, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi02 = LinearSpline(logr, xitmp);
    ComputeXiLM(1, 1, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi11 = LinearSpline(logr, xitmp);
    ComputeXiLM(2, 0, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi20 = LinearSpline(logr, xitmp);
    ComputeXiLM(2, 2, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi22 = LinearSpline(logr, xitmp);
    ComputeXiLM(3, 1, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi31 = LinearSpline(logr, xitmp);
    ComputeXiLM(4, 2, kmin, kmax, P, Nr, &r[0], &xitmp[0]);
    szapudi->xi42 = LinearSpline(logr, xitmp);
}

XiFunc SzapudiModel::GetXi() {
    SzapudiXi* szapudi = new SzapudiXi();
    szapudi->f = f;
    szapudi->partial_f = false;
    FillSzapudiXi(szapudi, 8192, 1e-2, 1e4, 0., 10., pk);
    return XiFunc(szapudi);
}

XiFunc SzapudiModel::GetXiDeriv(int n) {
    assert(0 <= n && n < Nbands+1);

    if(n < Nbands) {    // P(k) band parameter
        double kk[2] = { bands[n].min, bands[n].max };
        double pk[2] = { 1, 1 };
        Spline P = LinearSpline(2, kk, pk);

        SzapudiXi* szapudi = new SzapudiXi();
        szapudi->f = f;
        szapudi->partial_f = false;
        FillSzapudiXi(szapudi, 8192, 1e-2, 1e4, bands[n].min, bands[n].max, P);
        return XiFunc(szapudi);
    }
    else {              // anisotropy parameter
        SzapudiXi* szapudi = new SzapudiXi();
        szapudi->f = f;
        szapudi->partial_f = true;
        FillSzapudiXi(szapudi, 8192, 1e-2, 1e4, 0., 10., pk);
        return XiFunc(szapudi);
    }
}
