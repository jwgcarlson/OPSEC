/* TODO: - refactor ComputeSignalMatrixS() to avoid 6-deep nested for loop.
 *       - parallelize ComputeSignalMatrixS() in the same fashion as
 *         ComputeSignalMatrixC()
 *       - explain algorithm in more detail (in particular all the extra
 *         book-keeping involved to exploit azimuthal symmetry) */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
using std::vector;

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "sig.h"
#include "SeparationFunc.h"
#include "SelectionFunc.h"
#include "Survey.h"
#include "XiFunc.h"
#include "abn.h"
#include "cfg.h"
#include "cubature.h"
#include "opsec.h"
#include "rng.h"

#ifndef SIGNAL_MAXEVAL
#  define SIGNAL_MAXEVAL 50000000
#endif

template<class T> static inline T cb(T x) { return x*x*x; }
template<class T> static inline void swap(T& x, T& y) { T tmp = x; x = y; y = tmp; }

/* Parameters:
 * Ncells: size of global signal matrix
 * rows: array of row indices for the requested block (need not be consecutive)
 * cols: array of column indices for the requested block (need not be consecutive)
 * s: array of length at least nrows*lld
 * lld: leading dimension of local block s
 * xi: 2-point correlation function
 * ... : ... */
void ComputeSignalMatrixBlock(
        int Ncells, const vector<int>& rows, const vector<int>& cols,
        real* s, int lld,
        const XiFunc& xi, Survey* survey,
        const Cell* cells, int N1, int N2, int N3,
        double epsrel, double epsabs)
{
    if(lld <= 0)
        lld = (int) cols.size();
    else if(lld < (int) cols.size()) {
        opsec_error("ComputeSignalMatrixBlock: invalid array stride (%d,%zd,%zd)\n", lld, rows.size(), cols.size());
        opsec_abort(1);
    }
    int coordsys = survey->GetCoordinateSystem();
    if(coordsys == CoordSysCartesian)
        ComputeSignalMatrixBlockC(Ncells, (int) rows.size(), &rows[0], (int) cols.size(), &cols[0], s, lld, xi, survey, cells, N1, N2, N3, epsrel, epsabs);
    else if(coordsys == CoordSysSpherical)
        ComputeSignalMatrixBlockS(Ncells, (int) rows.size(), &rows[0], (int) cols.size(), &cols[0], s, lld, xi, survey, cells, N1, N2, N3, epsrel, epsabs);
}

/* Parameters:
 * n: size of global signal matrix
 * nrows: number of rows in the requested block
 * rows: indices of requested rows (might not be consecutive)
 * ncols: number of columns in the requested block
 * cols: indices of requested columns (might not be consecutive)
 * s: array of length at least nrows*lld
 * lld: leading dimension of local block s
 * xi: 2-point correlation function
 * ... : ... */
void ComputeSignalMatrixBlockC(
        int Ncells, int nrows, const int* rows, int ncols, const int* cols,
        real* s, int lld,
        const XiFunc& xi, Survey* survey,
        const Cell* cells, int Nx, int Ny, int Nz,
        double epsrel, double epsabs)
{
    /* The matrix element for cells a and b depends only on the relative
     * positions of the cells:
     *   kd = abs(da - db), ke = abs(ea - eb), kf = abs(fa - fb)
     * We keep a cache of already computed matrix elements and use it whenever
     * possible. */
    const real uninitialized = -1e100;  // unphysical value used to indicate an uninitialized element
    vector<real> cache(Nx*Ny*Nz, uninitialized);

    SeparationFunc sep = survey->GetSeparationFunction();

//    #pragma omp parallel for
    for(int i = 0; i < nrows; i++) {
        int a = rows[i];
        int da = cells[a].G / (Ny*Nz);
        int ea = (cells[a].G / Nz) % Ny;
        int fa = cells[a].G % Nz;
        for(int j = 0; j < ncols; j++) {
            int b = cols[j];
            int db = cells[b].G / (Ny*Nz);
            int eb = (cells[b].G / Nz) % Ny;
            int fb = cells[b].G % Nz;

            int kd = abs(da - db);
            int ke = abs(ea - eb);
            int kf = abs(fa - fb);
            int k = (kd*Ny + ke)*Nz + kf;

            real Q = cache[k];
            if(Q == uninitialized) {
                Q = (a <= b) ? ComputeSignalC(cells[a], cells[b], xi, sep, epsrel, epsabs)
                             : ComputeSignalC(cells[b], cells[a], xi, sep, epsrel, epsabs);
                cache[k] = Q;
            }

            s[i*lld + j] = sqrt(cells[a].Nbar * cells[b].Nbar) * Q;
        }
    }
}

void ComputeSignalMatrixBlockS(
        int Ncells, int nrows, const int* rows, int ncols, const int* cols,
        real* s, int lld,
        const XiFunc& xi, Survey* survey,
        const Cell* cells, int Nr, int Nmu, int Nphi,
        double epsrel, double epsabs)
{
    SeparationFunc sep = survey->GetSeparationFunction();

    const real uninitialized = -1e100;  // arbitrary magic number

    vector<int> rows_g2l(Ncells, -1), cols_g2l(Ncells, -1);
    for(int i = 0; i < nrows; i++)
        rows_g2l[rows[i]] = i;
    for(int j = 0; j < ncols; j++)
        cols_g2l[cols[j]] = j;

    /* Set local signal matrix values to uninitialized */
    for(int i = 0; i < nrows; i++)
        for(int j = 0; j < ncols; j++)
            s[i + j*lld] = uninitialized;

    /* Rearrange cells into a set of "arcs".  Each arc is a sequence of cells
     * at fixed r and mu, labeled by h = d*Nmu + e. */
    vector<vector<int> > arcs(Nr*Nmu);
    for(int a = 0; a < Ncells; a++) {
        int h = cells[a].G / Nphi;
        arcs[h].push_back(a);
    }

    vector<real> cache(Nphi);   // temporary list of stored values $Q_{ab}$

    /* Run ha over all arcs */
    for(int ha = 0; ha < Nr*Nmu; ha++) {
        int na = (int) arcs[ha].size();      // number of non-empty cells in arc ha
        if(na == 0)
            continue;
        int* acells = &arcs[ha][0];          // cells participating in arc ha

        /* Run hb over all arcs */
        for(int hb = 0; hb < Nr*Nmu; hb++) {
            int nb = (int) arcs[hb].size();  // number of non-empty cells in arc hb
            if(nb == 0)
                continue;
            int* bcells = &arcs[hb][0];      // cells participating in arc hb

            /* Clear list of temporary stored values */
            for(int k = 0; k < Nphi; k++)
                cache[k] = uninitialized;

            /* For fixed arcs ha and hb, iterate over cell pairs */
//            #pragma omp parallel for
            for(int ia = 0; ia < na; ia++) {
                int a = acells[ia];
                int i = rows_g2l[a];
                if(i < 0)
                    continue;

                for(int ib = 0; ib < nb; ib++) {
                    int b = bcells[ib];
                    int j = cols_g2l[b];
                    if(j < 0)
                        continue;

                    /* Compute matrix element only if it cannot be determined
                     * by symmetry with a previously computed element */
                    int fa = cells[a].G % Nphi;
                    int fb = cells[b].G % Nphi;
                    int k = abs(fa - fb);
                    real Q = cache[k];
                    if(Q == uninitialized) {
                        /* Update cache element */
                        Q = (a <= b) ? ComputeSignalS(cells[a], cells[b], xi, sep, epsrel, epsabs)
                                     : ComputeSignalS(cells[b], cells[a], xi, sep, epsrel, epsabs);
                        cache[k] = Q;
                    }
                    /* Set S_{ab} = \sqrt{\Nbar_a \Nbar_b} Q_{ab} */
                    s[i*lld + j] = sqrt(cells[a].Nbar * cells[b].Nbar) * Q;
                }
            }
        }
    }

    /* Sanity check: all local matrix elements should be set */
    for(int i = 0; i < nrows; i++)
        for(int j = 0; j < ncols; j++)
            assert(s[i*lld + j] != uninitialized);
}

#if 0
/* Compute the (unnormalized) signal matrix element $S_{ab}$, where a = c1.a
 * and b = c2.a, by performing a Monte Carlo integral over the 6-dimensional
 * volume $V_a \times V_b$.  The normalization factor $\sqrt{\Nbar_a \Nbar_b}$
 * is intentionally left off, so that the same element can be reused for other
 * pairs of cells that have different selection function values but the same
 * relative separation.  (Don't worry if this doesn't make sense.) */
double ComputeSignalS(const Cell& c1, const Cell& c2, const XiFunc& xi, const SeparationFunc& sep, double epsrel, double epsabs) {
    const int lg2n_min = 6, lg2n_max = 30;

    double rmin1cb = cb(c1.rmin), rmax1cb = cb(c1.rmax);
    double rmin2cb = cb(c2.rmin), rmax2cb = cb(c2.rmax);
    double phimin1 = c1.phimin, phimax1 = c1.phimax;
    double phimin2 = c2.phimin, phimax2 = c2.phimax;
    if(phimin1 > phimin2) {
        /* Swap phi boundaries if necessary (exploiting azimuthal symmetry of
         * cells), to ensure consistent numerical values for matrix elements
         * that should be identical. */
        swap(phimin1, phimin2);
        swap(phimax1, phimax2);
    }

    /* Initialize quasi-random sequence generator. */
    Sobol sobol;
    rng_sobol_init(&sobol, 6);
    double u[6];        // quasi-random point in the 6-D unit cube

    int n = 0;
    double fsum = 0.;
    double Q, oldQ = 0., dQ;
    Point p1, p2;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        while(n < (1 << lg2n)) {
            rng_sobol_get(&sobol, u);
            p1.r = cbrt((1-u[0]) * rmin1cb +  u[0] * rmax1cb);
            p2.r = cbrt((1-u[1]) * rmin2cb +  u[1] * rmax2cb);
            p1.mu =     (1-u[2]) * c1.mumin + u[2] * c1.mumax;
            p2.mu =     (1-u[3]) * c2.mumin + u[3] * c2.mumax;
            p1.phi =    (1-u[4]) * phimin1 +  u[4] * phimax1;
            p2.phi =    (1-u[5]) * phimin2 +  u[5] * phimax2;

            fsum += xi(p1, p2, sep);
            n++;
        }

        /* Declare convergence if the integral changes by only a small amount
         * after doubling the number of evaluation points. */
        Q = fsum/n;
        dQ = fabs(Q - oldQ);
        if(dQ < fmax(epsabs, fabs(Q)*epsrel))
            break;

        /* Otherwise continue with the next batch of evaluations. */
        oldQ = Q;
    }
    return Q;
}
#endif

real* ComputeSignalMatrixS(
    size_t n, size_t nloc, int amin, const Cell* cells,
    int Nr, int Nmu, int Nphi, const XiFunc& xi, Survey* survey,
    double epsrel, double epsabs)
{
    SeparationFunc sep = survey->GetSeparationFunction();

    const double uninitialized = -1e100;        // arbitrary magic number
    int amax = amin + nloc - 1;

    /* Allocate memory for local signal matrix */
    real* S = (real*) opsec_malloc(nloc*n*sizeof(real));
    if(S == NULL) {
        fprintf(stderr, "ComputeSignalMatrixS: could not allocate %zd bytes for local signal matrix\n", nloc*n*sizeof(real));
        return NULL;
    }
    opsec_info("Allocated %.1f MB for local signal matrix.\n", nloc*double(n*sizeof(real))/(1024.*1024.));

    /* Set local signal matrix values to uninitialized */
    for(ptrdiff_t i = 0; i < n*nloc; i++)
        S[i] = uninitialized;

    /* Rearrange cells into groups that share the same d and e */
    vector<vector<int> > decells(Nr*Nmu); // decells[h] is a list of all cells sharing the same d and e, where h = d*Nmu + e
    for(int a = 0; a < n; a++) {
        int h = cells[a].G / Nphi;
        decells[h].push_back(a);
    }

    /* Compute signal matrix $S_{ab}$, stored locally in row-major format.  Use
     * symmetries to avoid calculating elements that can be determined
     * otherwise.  */

    vector<double> cache(Nphi);   // temporary list of stored values $Q_{ab}$

    int ha, hb, a, b;
//    int da, ea, db, eb;
    int fa, fb, k;
    int na, nb, ia, ib;
    int *acells, *bcells;
    double Q;

    /* Iterate over (da,ea) pairs */
    for(ha = 0; ha < Nr*Nmu; ha++) {
//        da = ha / Nmu;
//        ea = ha % Nmu;
        na = (int) decells[ha].size();
        if(na == 0)
            continue;
        acells = &decells[ha][0];

        /* Iterate over (db,eb) pairs */
        for(hb = 0; hb < Nr*Nmu; hb++) {
//            db = hb / Nmu;
//            eb = hb % Nmu;
            nb = (int) decells[hb].size();
            if(nb == 0)
                continue;
            bcells = &decells[hb][0];

            /* Clear list of temporary stored values */
            for(k = 0; k < Nphi; k++)
                cache[k] = uninitialized;

            /* For fixed (da,ea) and (db,eb), iterate over cell pairs */
            for(ia = 0; ia < na; ia++) {
                a = acells[ia];
                if(a < amin || a > amax)
                    continue;

                for(ib = 0; ib < nb; ib++) {
                    b = bcells[ib];

                    /* Compute matrix element only if it cannot be determined
                     * by symmetry with a previously computed element */
                    if(amin <= b && b <= amax && S[a + (b-amin)*n] != uninitialized) {
                        /* Set S_{ab} equal to S_{ba} */
                        S[b + (a-amin)*n] = S[a + (b-amin)*n];
                    }
                    else {
                        fa = cells[a].G % Nphi;
                        fb = cells[b].G % Nphi;
                        k = abs(fa - fb);
                        if(cache[k] == uninitialized) {
                            /* Update cache element */
                            Q = (a <= b) ? ComputeSignalS(cells[a], cells[b], xi, sep, epsrel, epsabs)
                                         : ComputeSignalS(cells[b], cells[a], xi, sep, epsrel, epsabs);
                            cache[k] = Q;
                        }
                        /* Set S_{ab} = \sqrt{\Nbar_a \Nbar_b} Q_{ab} */
                        S[b + (a-amin)*n] = sqrt(cells[a].Nbar * cells[b].Nbar) * cache[k];
                    }
                }
            }
        }
    }

    /* Sanity check: all local matrix elements should be set */
    for(ptrdiff_t i = 0; i < n*nloc; i++) {
        assert(S[i] != uninitialized);
    }

    return S;
}

struct CubatureData {
    double min[6];
    double max[6];
    const XiFunc& xi;
    const SeparationFunc& sep;
};

static void cartesian_integrand(unsigned ndim, unsigned npt, const double* x, void* vdata, unsigned fdim, double* fval) {
#ifdef OPSEC_DEBUG
    assert(ndim == 6 && fdim == 1);
#endif
    CubatureData* data = (CubatureData*) vdata;
    double xmin1 = data->min[0], ymin1 = data->min[1], zmin1 = data->min[2],
           xmin2 = data->min[3], ymin2 = data->min[4], zmin2 = data->min[5];
    double xmax1 = data->max[0], ymax1 = data->max[1], zmax1 = data->max[2],
           xmax2 = data->max[3], ymax2 = data->max[4], zmax2 = data->max[5];
    const XiFunc& xi = data->xi;
    const SeparationFunc& sep = data->sep;
    #pragma omp parallel
    {
        Point p1, p2;
        #pragma omp for
        for(unsigned i = 0; i < npt; i++) {
            const double* u = &x[6*i];
            p1.x = (1 - u[0])*xmin1 + u[0]*xmax1;
            p1.y = (1 - u[1])*ymin1 + u[1]*ymax1;
            p1.z = (1 - u[2])*zmin1 + u[2]*zmax1;
            p2.x = (1 - u[3])*xmin2 + u[3]*xmax2;
            p2.y = (1 - u[4])*ymin2 + u[4]*ymax2;
            p2.z = (1 - u[5])*zmin2 + u[5]*zmax2;
            fval[i] = xi(p1, p2, sep);
        }
    }
}

double ComputeSignalC(const Cell& c1, const Cell& c2, const XiFunc& xi, const SeparationFunc& sep, double epsrel, double epsabs, double* err_) {
    CubatureData data = {
        { c1.xmin, c1.ymin, c1.zmin, c2.xmin, c2.ymin, c2.zmin },
        { c1.xmax, c1.ymax, c1.zmax, c2.xmax, c2.ymax, c2.zmax },
        xi,
        sep
    };
    double xmin[6] = { 0, 0, 0, 0, 0, 0 };
    double xmax[6] = { 1, 1, 1, 1, 1, 1 };
    unsigned maxeval = SIGNAL_MAXEVAL;
    double val, err;
    int okay = adapt_integrate_v(1, cartesian_integrand, (void*) &data,
                                 6, xmin, xmax,
                                 maxeval, epsrel, epsabs,
                                 &val, &err);
    if(err_ != NULL)
        *err_ = err;
    return val;
}

static void spherical_integrand(unsigned ndim, unsigned npt, const double* x, void* vdata, unsigned fdim, double* fval) {
#ifdef OPSEC_DEBUG
    assert(ndim == 6 && fdim == 1);
#endif
    CubatureData* data = (CubatureData*) vdata;
    const double rmin1cb = data->min[0], mumin1 = data->min[1], phimin1 = data->min[2],
           rmin2cb = data->min[3], mumin2 = data->min[4], phimin2 = data->min[5];
    const double rmax1cb = data->max[0], mumax1 = data->max[1], phimax1 = data->max[2],
           rmax2cb = data->max[3], mumax2 = data->max[4], phimax2 = data->max[5];
    const XiFunc& xi = data->xi;
    const SeparationFunc& sep = data->sep;
    #pragma omp parallel
    {
        Point p1, p2;
        #pragma omp for
        for(unsigned i = 0; i < npt; i++) {
            const double* u = &x[6*i];
            p1.r = cbrt((1-u[0]) * rmin1cb +  u[0] * rmax1cb);
            p2.r = cbrt((1-u[1]) * rmin2cb +  u[1] * rmax2cb);
            p1.mu =     (1-u[2]) * mumin1  +  u[2] * mumax1;
            p2.mu =     (1-u[3]) * mumin2  +  u[3] * mumax2;
            p1.phi =    (1-u[4]) * phimin1 +  u[4] * phimax1;
            p2.phi =    (1-u[5]) * phimin2 +  u[5] * phimax2;
            fval[i] = xi(p1, p2, sep);
        }
    }
}

double ComputeSignalS(const Cell& c1, const Cell& c2, const XiFunc& xi, const SeparationFunc& sep, double epsrel, double epsabs, double* err_) {
    CubatureData data = {
        { cb(c1.rmin), c1.mumin, c1.phimin, cb(c2.rmin), c2.mumin, c2.phimin },
        { cb(c1.rmax), c1.mumax, c1.phimax, cb(c2.rmax), c2.mumax, c2.phimax },
        xi,
        sep
    };
    double xmin[6] = { 0, 0, 0, 0, 0, 0 };
    double xmax[6] = { 1, 1, 1, 1, 1, 1 };
    unsigned maxeval = SIGNAL_MAXEVAL;
    double val, err;
    adapt_integrate_v(1, spherical_integrand, (void*) &data,
                      6, xmin, xmax,
                      maxeval, epsrel, epsabs,
                      &val, &err);
    if(err_ != NULL)
        *err_ = err;
    return val;
}

#if 0
double ComputeSignalC(const Cell& c1, const Cell& c2, const XiFunc& xi, const SeparationFunc& sep, double epsrel, double epsabs, int* neval) {
    const int lg2n_min = 3, lg2n_max = 25;

    /* Swap coordinate boundaries (exploiting translational symmetry) to
     * ensure consistent numerical values for matrix elements. */
    double xmin1 = c1.xmin, xmax1 = c1.xmax;
    double xmin2 = c2.xmin, xmax2 = c2.xmax;
    if(xmin1 > xmin2) {
        swap(xmin1, xmin2);
        swap(xmax1, xmax2);
    }
    double ymin1 = c1.ymin, ymax1 = c1.ymax;
    double ymin2 = c2.ymin, ymax2 = c2.ymax;
    if(ymin1 > ymin2) {
        swap(ymin1, ymin2);
        swap(ymax1, ymax2);
    }
    double zmin1 = c1.zmin, zmax1 = c1.zmax;
    double zmin2 = c2.zmin, zmax2 = c2.zmax;
    if(zmin1 > zmin2) {
        swap(zmin1, zmin2);
        swap(zmax1, zmax2);
    }

    /* Initialize quasi-random sequence generator. */
//    Sobol sobol;
//    rng_sobol_init(&sobol, 6);
    rng_state rng;
    rng_init_r(&rng, 72);
    double u[6];        // quasi-random point in the 6-D unit cube

    int n = 0;
    double fsum = 0.;
    double Q, oldQ = 0.;
    Point p1, p2;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        int ngoal = (1 << lg2n);
        while(n < ngoal) {
//            rng_sobol_get(&sobol, u);
            for(int i = 0; i < 6; i++)
                u[i] = rng_uniform_r(&rng);
            p1.x = (1 - u[0])*xmin1 + u[0]*xmax1;
            p1.y = (1 - u[1])*ymin1 + u[1]*ymax1;
            p1.z = (1 - u[2])*zmin1 + u[2]*zmax1;
            p2.x = (1 - u[3])*xmin2 + u[3]*xmax2;
            p2.y = (1 - u[4])*ymin2 + u[4]*ymax2;
            p2.z = (1 - u[5])*zmin2 + u[5]*zmax2;

            fsum += xi(p1, p2, sep);
            n++;
        }

        /* Declare convergence if the integral changes by only a small amount
         * after doubling the number of evaluation points. */
        Q = fsum/n;
        double dQ = fabs(Q - oldQ);
        if(dQ < fmax(epsabs, fabs(Q)*epsrel))
            break;

        /* Otherwise continue with the next batch of evaluations. */
        oldQ = Q;

        if(lg2n == lg2n_max)
            fprintf(stderr, "ComputeSignalC: max evaluations reached\n");
    }
    if(neval)
        *neval = n;
    return Q;
}
#endif

real* ComputeSignalMatrixC(
    size_t n, size_t nloc, int amin,
    const Cell* cells,
    int Nx, int Ny, int Nz,
    const XiFunc& xi, Survey* survey,
    double epsrel, double epsabs)
{
    SeparationFunc sep = survey->GetSeparationFunction();

    /* Allocate memory for local signal matrix */
    real* S = (real*) opsec_malloc(nloc*n*sizeof(real));
    if(S == NULL) {
        fprintf(stderr, "ComputeSignalMatrixC: could not allocate %zd bytes for local signal matrix\n", nloc*n*sizeof(real));
        return NULL;
    }

    /* The matrix element for cells a and b depends only on the relative
     * positions of the cells:
     *   Dd = abs(da - db), De = abs(ea - eb), Df = abs(fa - fb)
     * We keep a cache of already computed matrix elements and use it whenever
     * possible. */
    const real uninitialized = -1e100;  // unphysical value used to indicate an uninitialized element
    vector<real> cache(Nx*Ny*Nz, uninitialized);

    int aloc, a, da, ea, fa;
    int       b, db, eb, fb;
    int Dd, De, Df;
    real Q;
    #pragma omp parallel for private(aloc,a,da,ea,fa,b,db,eb,fb,Dd,De,Df,Q)
    for(aloc = 0; aloc < nloc; aloc++) {
        a = aloc + amin;
        da = cells[a].G / (Ny*Nz);
        ea = (cells[a].G / Nz) % Ny;
        fa = cells[a].G % Nz;
        for(b = 0; b < n; b++) {
            db = cells[b].G / (Ny*Nz);
            eb = (cells[b].G / Nz) % Ny;
            fb = cells[b].G % Nz;

            Dd = abs(da - db);
            De = abs(ea - eb);
            Df = abs(fa - fb);

            Q = cache[(Dd*Ny + De)*Nz + Df];
            if(Q == uninitialized) {
                Q = (a <= b) ? ComputeSignalC(cells[a], cells[b], xi, sep, epsrel, epsabs)
                             : ComputeSignalC(cells[b], cells[a], xi, sep, epsrel, epsabs);
                cache[(Dd*Ny + De)*Nz + Df] = Q;
                /* Even though cache is shared across threads, and race
                 * conditions are likely to occur, the code above should be
                 * safe.  In the worst case, the same matrix element will be
                 * computed by more than one thread.  This is inefficient, but
                 * it shouldn't happen very often, and in any case it will not
                 * change the final answer. */
            }

            S[aloc*n + b] = sqrt(cells[a].Nbar * cells[b].Nbar) * Q;
        }
    }

    return S;
}

FILE* ReadModesHeader(const char* modefile, int* Nmodes_, int* Ncells_) {
    /* Read in spectrum from file */
    FILE* fp = fopen(modefile, "r");
    if(fp == NULL) {
        fprintf(stderr, "ReadModesHeader: could not open file '%s'\n", modefile);
        return NULL;
    }

    size_t n, size;
    Config opts = cfg_new();
    int err = abn_read_header(fp, &n, &size, NULL, NULL, opts);
    int Nmodes = cfg_get_int(opts, "Nmodes");
    int Ncells = cfg_get_int(opts, "Ncells");
    cfg_destroy(opts);

    if(err || Nmodes <= 0 || Ncells <= 0 || n != Nmodes*Ncells || size != sizeof(double)) {
        fprintf(stderr, "ReadModesHeader: error reading '%s'\n", modefile);
        fclose(fp);
        return NULL;
    }

    if(Nmodes_) *Nmodes_ = Nmodes;
    if(Ncells_) *Ncells_ = Ncells;
    return fp;
}

real* ReadModes(const char* modefile, int* Nmodes_, int* Ncells_) {
    int Nmodes = -1;
    int Ncells = -1;
    real* modes = NULL;
    int nread;

    int me = 0, usempi = 0;
#ifdef OPSEC_USE_MPI
    MPI_Initialized(&usempi);
    if(usempi)
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    /* Read in spectrum from file on the root process */
    if(me == 0) {
        FILE* fp = ReadModesHeader(modefile, &Nmodes, &Ncells);
        if(fp != NULL) {
            modes = (real*) opsec_malloc(Nmodes*Ncells*sizeof(real));
            if(modes == NULL)
                fprintf(stderr, "ReadModes: could not allocate %zd bytes of memory\n", Nmodes*Ncells*sizeof(real));
            else {
                if(fread(modes, sizeof(real), Nmodes*Ncells, fp) != Nmodes*Ncells) {
                    fprintf(stderr, "ReadModes: error reading data from '%s'\n", modefile);
                    free(modes);
                    modes = NULL;
                }
            }
            fclose(fp);
        }
    }

#ifdef OPSEC_USE_MPI
    /* Broadcast spectrum to other processes */
    if(usempi) {
        MPI_Bcast(&Nmodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Ncells, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(Nmodes > 0) {
            if(me != 0)
                modes = (real*) opsec_malloc(Nmodes*Ncells*sizeof(real));
            MPI_Bcast(modes, Nmodes*Ncells, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
        }
    }
#endif

    if(Nmodes_) *Nmodes_ = Nmodes;
    if(Ncells_) *Ncells_ = Ncells;
    return modes;
}
