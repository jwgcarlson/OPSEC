#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "BoxSurvey.h"
#include "SelectionFunc.h"
#include "SeparationFunc.h"
#include "Spline.h"
#include "abn.h"


/* Separation between two points in a periodic box. */
struct BoxSeparationFuncImpl : public SeparationFuncImpl {
    double L;   // box size
    Point O;    // origin

    BoxSeparationFuncImpl(double L_, const Point& O_) : L(L_), O(O_) {}

    /* Separation between two points x1 and x2 on the periodic interval [0,L] */
    double linsep(double x1, double x2) {
        double dx = x2 - x1;
        if(dx > L/2)
            dx -= L;
        else if(dx < -L/2)
            dx += L;
        return dx;
    }

    double r(const Point& p1, const Point& p2) {
        double rx = linsep(p1.x, p2.x);
        double ry = linsep(p1.y, p2.y);
        double rz = linsep(p1.z, p2.z);
        return sqrt(rx*rx + ry*ry + rz*rz);
    }

    void rmu(const Point& p1, const Point& p2, double& r, double& mu) {
        double rx = linsep(p1.x, p2.x);
        double ry = linsep(p1.y, p2.y);
        double rz = linsep(p1.z, p2.z);
        r = sqrt(rx*rx + ry*ry + rz*rz);
        mu = rz/r;
    }

    void rab(const Point& p1, const Point& p2, double& r, double& a, double& b) {
        double x1 = p1.x - O.x;
        double y1 = p1.y - O.y;
        double z1 = p1.z - O.z;
        double x2 = p2.x - O.x;
        double y2 = p2.y - O.y;
        double z2 = p2.z - O.z;
        a = sqrt(x1*x1 + y1*y1 + z1*z1);
        b = sqrt(x2*x2 + y2*y2 + z2*z2);
        r = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
    }
};

/* Constant selection function. */
struct ConstantSelectionFuncImpl : public SelectionFuncImpl {
    double const_nbar;          // constant value for $\bar{n}(\vec{r})$

    ConstantSelectionFuncImpl(double const_nbar_) {
        const_nbar = const_nbar_;
        assert(const_nbar > 0);
    }

    ~ConstantSelectionFuncImpl() {}

    double evaluate(double, double, double) {
        return const_nbar;
    }
};


BoxSurvey::BoxSurvey(Config cfg) {
    galfile = cfg_get(cfg, "galfile");
    nbar = cfg_get_double(cfg, "nbar");
    L = cfg_get_double(cfg, "L");
}

BoxSurvey::~BoxSurvey() {
}

int BoxSurvey::GetCoordinateSystem() const {
    return CoordSysCartesian;
}

SeparationFunc BoxSurvey::GetSeparationFunction() {
    Point origin = { L/2, L/2, L/2 };
    return SeparationFunc(new BoxSeparationFuncImpl(L, origin));
}

SelectionFunc BoxSurvey::GetSelectionFunction() {
    return new ConstantSelectionFuncImpl(nbar);
}

void BoxSurvey::GetGalaxies(std::vector<Galaxy>& gals) {
    size_t norig = gals.size();

    FILE* fgals = fopen(galfile.c_str(), "r");
    if(fgals == NULL) {
        fprintf(stderr, "BoxSurvey: could not open file '%s'\n", galfile.c_str());
        return;
    }

    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    if(abn_read_header(fgals, &n, &size, &endian, fmt, NULL) != 0
       || size != sizeof(Galaxy)
       || strcmp(fmt, "4f") != 0)
    {
        fprintf(stderr, "BoxSurvey: error reading galaxies from '%s'\n", galfile.c_str());
        fclose(fgals);
        return;
    }

    /* Make room for n additional Galaxy objects */
    gals.resize(norig + n);

    size_t nread = fread((void*) &gals[norig], size, n, fgals);
    if(nread != n)
        fprintf(stderr, "BoxSurvey: expecting %zd galaxies from '%s', got %zd\n", n, galfile.c_str(), nread);

    fclose(fgals);
}
