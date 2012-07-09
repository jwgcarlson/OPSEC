#include <cmath>
#include <cstring>

#include "Survey.h"
#include "opsec.h"

/* Return a representative of the periodic coordinate x in the interval [0,L). */
static inline double periodic(double x, double L = 2*M_PI) {
    return x - floor(x/L) * L;
}

Survey::Survey(int coordsys, Config cfg)
    : coordsys(coordsys)
{
    if(cfg) {
        /* Read common survey parameters */
        if(coordsys == CoordSysCartesian) {
            if(char* missing = cfg_missing_keys(cfg, "XMin,XMax,YMin,YMax,ZMin,ZMax,Nx,Ny,Nz")) {
                opsec_error("Survey: missing configuration options %s\n", missing);
                opsec_abort(1);
            }
            XMin = cfg_get_double(cfg, "XMin");
            XMax = cfg_get_double(cfg, "XMax");
            YMin = cfg_get_double(cfg, "YMin");
            YMax = cfg_get_double(cfg, "YMax");
            ZMin = cfg_get_double(cfg, "ZMin");
            ZMax = cfg_get_double(cfg, "ZMax");
            Nx = cfg_get_int(cfg, "Nx");
            Ny = cfg_get_int(cfg, "Ny");
            Nz = cfg_get_int(cfg, "Nz");
        }
        else if(coordsys == CoordSysSpherical) {
            if(char* missing = cfg_missing_keys(cfg, "RMin,RMax,MuMin,MuMax,PhiMin,PhiMax,Nr,Nmu,Nphi")) {
                opsec_error("Survey: missing configuration options %s\n", missing);
                opsec_abort(1);
            }
            RMin = cfg_get_double(cfg, "RMin");
            RMax = cfg_get_double(cfg, "RMax");
            MuMin = cfg_get_double(cfg, "MuMuin");
            MuMax = cfg_get_double(cfg, "MuMax");
            PhiMin = cfg_get_double(cfg, "PhiMin");
            PhiMax = cfg_get_double(cfg, "PhiMax");
            Nr = cfg_get_int(cfg, "Nr");
            Nmu = cfg_get_int(cfg, "Nmu");
            Nphi = cfg_get_int(cfg, "Nphi");
        }
    }
    else {
        coordsys = CoordSysUndefined;
        Min1 = Min2 = Min3 = 0;
        Max1 = Max2 = Max3 = 1;
        N1 = N2 = N3 = 1;
    }

    /* Validate parameters */
    opsec_assert(Min1 < Max1);
    opsec_assert(Min2 < Max2);
    opsec_assert(Min3 < Max3);
    opsec_assert(N1 > 0);
    opsec_assert(N2 > 0);
    opsec_assert(N3 > 0);
}

Survey::~Survey() {
}

int Survey::GetCoordinateSystem() const {
    return coordsys;
}

Cell Survey::CreateEmptyCell(int d1, int d2, int d3) const {
    Cell c;
    c.a = -1;
    c.d1 = d1;
    c.d2 = d2;
    c.d3 = d3;
    c.min1 = Min1 +     d1*(Max1 - Min1)/N1;
    c.max1 = Min1 + (d1+1)*(Max1 - Min1)/N1;
    c.min2 = Min2 +     d2*(Max2 - Min2)/N2;
    c.max2 = Min2 + (d2+1)*(Max2 - Min2)/N2;
    c.min3 = Min3 +     d3*(Max3 - Min3)/N3;
    c.max3 = Min3 + (d3+1)*(Max3 - Min3)/N3;
    if(coordsys == CoordSysCartesian)
        c.Veff = (c.xmax - c.xmin) * (c.ymax - c.ymin) * (c.zmax - c.zmin);
    else if(coordsys == CoordSysSpherical)
        c.Veff = (pow(c.rmax, 3) - pow(c.rmin, 3))/3. * (c.mumax - c.mumin) * (c.phimax - c.phimin);
    else
        c.Veff = 0;
    c.Nbar = 0;
    return c;
}

int Survey::GetGridIndex(double x1, double x2, double x3) const {
    int d1 = (int) floor(N1 * (x1 - Min1)/(Max1 - Min1));
    int d2 = (int) floor(N2 * (x2 - Min2)/(Max2 - Min2));
    int d3 = (int) floor(N3 * (x3 - Min3)/(Max3 - Min3));
    if(d1 < 0 || d1 >= N1 || d2 < 0 || d2 >= N2 || d3 < 0 || d3 >= N3)
        return -1;
    else
        return (d1*N2 + d2)*N3 + d3;
}

Config Survey::GetConfigurationOptions() const {
    Config opts = cfg_new();
    if(coordsys == CoordSysCartesian) {
        cfg_set(opts, "coordsys", "cartesian");
        cfg_set_double(opts, "XMin", XMin);
        cfg_set_double(opts, "XMax", XMax);
        cfg_set_double(opts, "YMin", YMin);
        cfg_set_double(opts, "YMax", YMax);
        cfg_set_double(opts, "ZMin", ZMin);
        cfg_set_double(opts, "ZMax", ZMax);
        cfg_set_int(opts, "Nx", Nx);
        cfg_set_int(opts, "Ny", Ny);
        cfg_set_int(opts, "Nz", Nz);
    }
    else if(coordsys == CoordSysSpherical) {
        cfg_set(opts, "coordsys", "spherical");
        cfg_set_double(opts, "RMin", RMin);
        cfg_set_double(opts, "RMax", RMax);
        cfg_set_double(opts, "MuMin", MuMin);
        cfg_set_double(opts, "MuMax", MuMax);
        cfg_set_double(opts, "PhiMin", PhiMin);
        cfg_set_double(opts, "PhiMax", PhiMax);
        cfg_set_int(opts, "Nr", Nr);
        cfg_set_int(opts, "Nmu", Nmu);
        cfg_set_int(opts, "Nphi", Nphi);
    }
    return opts;
}


/* If you create your own Survey sub-class, put the necessary #include
 * statement here and add a corresponding HANDLE_SURVEY() call below.  Your
 * sub-class's constructor must take a Config object as its only argument. */

#include "BoxSurvey.h"
#include "SphericalSurvey.h"

Survey* InitializeSurvey(Config cfg) {
    if(!cfg_has_key(cfg, "survey.type")) {
        opsec_error("InitializeSurvey: missing survey.type configuration option\n");
        return NULL;
    }

    const char* type = cfg_get(cfg, "survey.type");
    Config surveycfg = cfg_new_sub(cfg, "survey.");
    Survey* s = NULL;

#define HANDLE_SURVEY(NAME) \
    if(strcmp(type, #NAME) == 0) { \
        s = new NAME(surveycfg); \
    }

    HANDLE_SURVEY(BoxSurvey)
    HANDLE_SURVEY(SphericalSurvey)
#undef HANDLE_SURVEY

    if(s == NULL) {
        opsec_error("InitializeSurvey: unrecognized survey type '%s'\n", type);
    }

    cfg_destroy(surveycfg);
    return s;
}
