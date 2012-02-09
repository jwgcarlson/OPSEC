/* TODO: - refactor ComputeSignalMatrixS() to avoid 6-deep nested for loop.
 *       - parallelize ComputeSignalMatrixS() in the same fashion as
 *         ComputeSignalMatrixC()
 *       - explain algorithm in more detail (in particular all the extra
 *         book-keeping involved to exploit azimuthal symmetry) */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
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

#include "abn.h"
#include "cfg.h"
#include "rng.h"
#include "sig.h"

template<typename T> static inline T sq(T x) { return x*x; }
template<typename T> static inline T cb(T x) { return x*x*x; }
template<typename T> static inline void swap(T& x, T& y) { T tmp = x; x = y; y = tmp; }


/* Compute the (unnormalized) signal matrix element $S_{ab}$, where a = c1.a
 * and b = c2.a, by performing a Monte Carlo integral over the 6-dimensional
 * volume $V_a \times V_b$.  The normalization factor $\sqrt{\Nbar_a \Nbar_b}$
 * is intentionally left off, so that the same element can be reused for other
 * pairs of cells that have different selection function values but the same
 * relative separation.  (Don't worry if this doesn't make sense.) */
double ComputeSignalS(const Cell& c1, const Cell& c2, const XiFunc& xi, double epsrel, double epsabs) {
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
    double p[6];
    rng_sobol_get(&sobol, p);   // throw away first point, which always gives r = 0

    int n = 0;
    double fsum = 0.;
    double Q, oldQ = 0., dQ;
    double r1, r2, mu1, mu2, phi1, phi2, r, cosgamma;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        while(n < (1 << lg2n)) {
            rng_sobol_get(&sobol, p);

            r1 = cbrt((1-p[0]) * rmin1cb +  p[0] * rmax1cb);
            r2 = cbrt((1-p[1]) * rmin2cb +  p[1] * rmax2cb);
            mu1 =     (1-p[2]) * c1.mumin + p[2] * c1.mumax;
            mu2 =     (1-p[3]) * c2.mumin + p[3] * c2.mumax;
            phi1 =    (1-p[4]) * phimin1 +  p[4] * phimax1;
            phi2 =    (1-p[5]) * phimin2 +  p[5] * phimax2;
            cosgamma = sqrt((1 - mu1*mu1)*(1 - mu2*mu2)) * cos(phi1 - phi2) + mu1*mu2;
            r = sqrt(fmax(r1*r1 + r2*r2 - 2*r1*r2*cosgamma, 0.0));

            fsum += xi(r, r1, r2);
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

real* ComputeSignalMatrixS(
    size_t n, size_t nloc, int amin, const Cell* cells,
    int Nr, int Nmu, int Nphi, const XiFunc& xi,
    double epsrel, double epsabs)
{
    const double uninitialized = -1e100;        // arbitrary magic number
    int amax = amin + nloc - 1;

    /* Allocate memory for local signal matrix */
    real* S = (real*) opsec_malloc(nloc*n*sizeof(real));
    if(S == NULL) {
        fprintf(stderr, "ComputeSignalMatrixS: could not allocate %zd bytes for local signal matrix\n", nloc*n*sizeof(real));
        return NULL;
    }
    printf("Allocated %.1f MB for local signal matrix.\n", nloc*double(n*sizeof(real))/(1024.*1024.));
    fflush(stdout);

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
            #pragma omp parallel for private(ia,ib,a,b,fa,fb,k,Q)
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
                            Q = (a <= b) ? ComputeSignalS(cells[a], cells[b], xi, epsrel, epsabs)
                                         : ComputeSignalS(cells[b], cells[a], xi, epsrel, epsabs);
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

static inline double norm(double dx, double dy, double dz) {
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double ComputeSignalC(const Cell& c1, const Cell& c2, const XiFunc& xi, double epsrel, double epsabs) {
    const int lg2n_min = 6, lg2n_max = 30;

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
    Sobol sobol;
    rng_sobol_init(&sobol, 6);
    double p[6];
    rng_sobol_get(&sobol, p);   // throw away first point, which always gives r = 0

    int n = 0;
    double fsum = 0.;
    double Q, oldQ = 0., dQ;
    double x1, y1, z1, x2, y2, z2, r, r1, r2;
    for(int lg2n = lg2n_min; lg2n <= lg2n_max; lg2n++) {
        while(n < (1 << lg2n)) {
            rng_sobol_get(&sobol, p);

            x1 = (1 - p[0])*xmin1 + p[0]*xmax1;
            y1 = (1 - p[1])*ymin1 + p[1]*ymax1;
            z1 = (1 - p[2])*zmin1 + p[2]*zmax1;
            x2 = (1 - p[3])*xmin2 + p[3]*xmax2;
            y2 = (1 - p[4])*ymin2 + p[4]*ymax2;
            z2 = (1 - p[5])*zmin2 + p[5]*zmax2;
            r1 = norm(x1, y1, z1);
            r2 = norm(x2, y2, z2);
            r = norm(x2 - x1, y2 - y1, z2 - z1);

            fsum += xi(r, r1, r2);
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

real* ComputeSignalMatrixC(
    size_t n, size_t nloc, int amin,
    const Cell* cells,
    int Nx, int Ny, int Nz,
    const XiFunc& xi,
    double epsrel, double epsabs)
{
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
                Q = (a <= b) ? ComputeSignalC(cells[a], cells[b], xi, epsrel, epsabs)
                             : ComputeSignalC(cells[b], cells[a], xi, epsrel, epsabs);
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

#if 0
real* ReadModes(const char* modefile, int* Nmodes_, int* Ncells_) {
    int Nmodes = -1;
    int Ncells = -1;
    real* modes = NULL;

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

    if(Nmodes_) *Nmodes_ = Nmodes;
    if(Ncells_) *Ncells_ = Ncells;
    return modes;
}
#endif

real* ReadModes(const char* modefile, int* Nmodes_, int* Ncells_) {
    int Nmodes = -1;
    int Ncells = -1;
    real* modes = NULL;
    int nread;

    int nprocs = 1, me = 0, usempi = 0;
#ifdef OPSEC_USE_MPI
    MPI_Initialized(&usempi);
    if(usempi) {
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
    }
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
