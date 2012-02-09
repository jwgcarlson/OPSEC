/* gmock
 *
 * Generate a mock catalog of "galaxies" by taking a Poisson sampling of a
 * Gaussian random field.
 *
 * Details:
 * - We start with a periodic cubical volume [0,L]^3, partitioned into an NxNxN
 *   grid.  We use the input power spectrum to generate a Gaussian random field
 *   $\delta(\vec{x})$ on the grid, by randomly drawing complex Fourier
 *   coefficients $\tilde{\delta}_{\vec{k}}$ with mean zero and variance
 *   $P(k)$, then performing a FFT.  (More specifically, the real and imaginary
 *   parts of $\tilde{\delta}_{\vec{k}}$ are drawn independently with variance
 *   $P(k)/2$, and reality conditions are enforced.)
 * - The density contrast is made into a number density field,
 *     n(\vec{x}) = \bar{n}_0 [1 + b \delta(\vec{x})]
 *   where $\bar{n}_0$ is the mean galaxy density and $b$ is a linear bias
 *   factor.  At grid points where this results in a negative number density,
 *   $n(\vec{x})$ is set to zero.
 * - Within each grid cell, galaxies are thrown down according to a Poisson
 *   point process.  That is, if the mean number density within cell $i$ (i.e.
 *   the average of $n(\vec{x})$ at all 8 corners) is $\bar{n}_i$, and the
 *   volume of the cell is $dV = (L/N)^3$, then the number of galaxies $N_i$
 *   to lay down within the cell is a Poisson random variable with mean
 *   $\lambda_i = \bar{n}_i dV$.  The positions of these $N_i$ galaxies are
 *   chosen in accordance with the number density within the grid cell, which
 *   is treated as a trilinear interpolation of the values $n(\vec{x})$ at the
 *   corners.
 * - Redshift distortions are performed if the anisotropy parameter f is
 *   nonzero.  In this case there are two options: either the observer is
 *   placed at a finite position (x,y,z), or the observer is placed at infinity
 *   with the line-of-sight direction along (nx,ny,nz).  The config parameter
 *   'observer' is used to specify the observer position.  Three possible
 *   choices, and their interpretations, are given below:
 *      observer = 0 0 -1 0             # line-of-sight along +z direction
 *      observer = 500 500 500 1        # observer placed at (x,y,z) = (500,500,500)
 *      observer = 500 500 500          # same (w defaults to 1)
 *   In other words, the observer position is specified in homogeneous
 *   coordinates (x,y,z,w), with the parameter w indicating either a finite
 *   observer position (w = 1, the default) or an observer at infinity (w = 0).
 *   (Technically, if w is nonzero, then the observer is placed at (x/w,y/w,z/w)).
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <vector>
using std::vector;

#include "abn.h"
#include "array.h"
#include "cfg.h"
#include "particle.h"
#include "rng.h"
#include "sampling.h"
#include "Spline.h"

struct Galaxy {
    float x, y, z;      // galaxy position (in real space or redshift space)
    float w;            // galaxy weight (always 1 for Gaussian mock galaxies)
};

/* Information describing a computational slab (i.e. the portion of the
 * computational box that belongs to a single process) */
struct Slab {
    double L;                   // dimensions of computational box
    int N;                      // number of grid cells in computational box (in each direction)
    int nxloc, ny, nz, ixmin;   // number of grid cells in slab, and position of slab within grid
    Grid rho, vx, vy, vz;       // density and velocity fields
};

const unsigned int DEFAULT_SEED = 1776;

/* Inverse CIC: return the value of the field at the point (x,y,z) \in [0,1]^3 */
static inline double icic(double a000, double a001, double a010, double a011,
                          double a100, double a101, double a110, double a111,
                          double x, double y, double z)
{
    return a000*(1-x)*(1-y)*(1-z) + a001*(1-x)*(1-y)*z + a010*(1-x)*y*(1-z) + a011*(1-x)*y*z
         + a100*x*(1-y)*(1-z)     + a101*x*(1-y)*z     + a110*x*y*(1-z)     + a111*x*y*z;
}

int sample_galaxies(Slab& slab, double* obs, vector<Galaxy>& galaxies) {
    double L = slab.L;
    int N = slab.N;
    int nxloc = slab.nxloc, ny = slab.ny, nz = slab.nz, ixmin = slab.ixmin;
    double dL = L/N;            // grid cell size
    double dV = dL*dL*dL;       // grid cell volume
    Grid& rho = slab.rho;
    Grid& vx = slab.vx;
    Grid& vy = slab.vy;
    Grid& vz = slab.vz;

    /* Whether or not to put the observer at infinity */
    bool distant_observer = (obs[3] == 0.);
    /* Line-of-sight direction in the distant observer case */
    double dx = -obs[0], dy = -obs[1], dz = -obs[2];
    /* Observer position otherwise */
    double ox = obs[0]/obs[3], oy = obs[1]/obs[3], oz = obs[2]/obs[3];

    int ntotal = 0;

    /* Iterate over grid cells within this slab */
    double x, y, z;
    double gvx, gvy, gvz;
    double gvn, norm;
    for(int ixloc = 0; ixloc < nxloc; ixloc++) {
        int ix = ixmin + ixloc;
        for(int iy = 0; iy < ny; iy++) {
            for(int iz = 0; iz < nz; iz++) {
                /* Compute mean density over the grid cell */
                double a000 = rho(ixloc,iy,iz);
                double a001 = rho(ixloc,iy,(iz+1)%nz);
                double a010 = rho(ixloc,(iy+1)%ny,iz);
                double a011 = rho(ixloc,(iy+1)%ny,(iz+1)%nz);
                double a100 = rho(ixloc+1,iy,iz);
                double a101 = rho(ixloc+1,iy,(iz+1)%nz);
                double a110 = rho(ixloc+1,(iy+1)%ny,iz);
                double a111 = rho(ixloc+1,(iy+1)%ny,(iz+1)%nz);
                double rhobar = (a000 + a001 + a010 + a011 + a100 + a101 + a110 + a111)/8;

                /* Decide how many new galaxies to add to this volume element */
                int numnew = (int)rng_poisson(rhobar*dV);
                double fx, fy, fz;      // position of particle within grid cell, (fx,fy,fz) \in [0,1]^3
                Galaxy g;
                g.w = 1;
                for(int i = 0; i < numnew; i++) {
                    /* Choose where within the grid cell to place the particle */
                    sample3d(a000, a001, a010, a011, a100, a101, a110, a111, fx, fy, fz);
                    x = (ix + fx)*dL;
                    y = (iy + fy)*dL;
                    z = (iz + fz)*dL;

                    /* Compute the velocity at that position */
                    gvx = icic(vx(ixloc,iy,iz), vx(ixloc,iy,(iz+1)%N), vx(ixloc,(iy+1)%N,iz), vx(ixloc,(iy+1)%N,(iz+1)%N), vx(ixloc+1,iy,iz), vx(ixloc+1,iy,(iz+1)%N), vx(ixloc+1,(iy+1)%N,iz), vx(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);
                    gvy = icic(vy(ixloc,iy,iz), vy(ixloc,iy,(iz+1)%N), vy(ixloc,(iy+1)%N,iz), vy(ixloc,(iy+1)%N,(iz+1)%N), vy(ixloc+1,iy,iz), vy(ixloc+1,iy,(iz+1)%N), vy(ixloc+1,(iy+1)%N,iz), vy(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);
                    gvz = icic(vz(ixloc,iy,iz), vz(ixloc,iy,(iz+1)%N), vz(ixloc,(iy+1)%N,iz), vz(ixloc,(iy+1)%N,(iz+1)%N), vz(ixloc+1,iy,iz), vz(ixloc+1,iy,(iz+1)%N), vz(ixloc+1,(iy+1)%N,iz), vz(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);

                    /* Compute redshift-space position */
                    if(distant_observer) {
                        g.x = x + dx*gvx;
                        g.y = y + dy*gvy;
                        g.z = z + dz*gvz;
                    }
                    else {
                        dx = x - ox;
                        dy = y - oy;
                        dz = z - oz;
                        norm = sqrt(dx*dx + dy*dy + dz*dz);
                        dx /= norm;
                        dy /= norm;
                        dz /= norm;
                        gvn = dx*gvx + dy*gvy + dz*gvz;
                        g.x = x + gvn*dx;
                        g.y = y + gvn*dy;
                        g.z = z + gvn*dz;
                    }

                    /* Add the particle */
                    galaxies.push_back(g);
                    ntotal++;
                }
            }
        }
    }

    return ntotal;
}

/* CIC windowing function, with $\vec{k} = (2\pi/N) \vec{n}$ given in grid
 * coordinates. */
static inline double W(double kx, double ky, double kz) {
    double wx = (kx == 0) ? 1 : sin(kx/2)/(kx/2);
    double wy = (ky == 0) ? 1 : sin(ky/2)/(ky/2);
    double wz = (kz == 0) ? 1 : sin(kz/2)/(kz/2);
    return pow(wx*wy*wz, 2);
}

#ifdef HAVE_MPI
void ComputePowerSpectrum(MPI_Comm comm, vector<Galaxy> galaxies, Slab& slab, const char* filename) {
    int me, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    double L = slab.L;
    int N = slab.N;
    int nxloc = slab.nxloc;
    int ixmin = slab.ixmin;
    Grid& rho = slab.rho;

    FILE* fout;
    if(me == 0) {
        printf("Writing measured P(k) to '%s'\n", filename);
        fout = fopen(filename, "w");
    }

    /* Define binning scheme */
    int nbins = 20;
    double kmin = 0.0;
    double kmax = 0.2;
    array power = array::zeros(nbins);      // sum of $|\delta_{\vec{k}}|^2$ within bin
    array count = array::zeros(nbins);      // number of modes contributing to bin

    /* Clear density field */
    rho.zero();

    /* Lay down galaxies on the grid */
    double x, y, z;
    int ixloc, ix, iy, iz;
    double fx, fy, fz;
    int ngals = (int)galaxies.size();
    for(int i = 0; i < ngals; i++) {
        /* Transform physical coordinates to grid coordinates (0 <= x < N), using periodicity */
        x = fmod(N*(galaxies[i].x/L + 1), N);
        y = fmod(N*(galaxies[i].y/L + 1), N);
        z = fmod(N*(galaxies[i].z/L + 1), N);

        /* Assign galaxy to grid using CIC */
        ix = (int)x;
        iy = (int)y;
        iz = (int)z;
        fx = x - ix;
        fy = y - iy;
        fz = z - iz;
        ixloc = ix - ixmin;
        if(!(0 <= ixloc && ixloc < N)) {
            printf("ix = %d, ixloc = %d, N = %d, me = %d\n", ix, ixloc, N, me);
        }
        assert(0 <= ixloc && ixloc < N);        // all galaxies for this process should belong to its own slab
        rho(ixloc  , iy      , iz      ) += (1-fx)*(1-fy)*(1-fz);
        rho(ixloc  , iy      , (iz+1)%N) += (1-fx)*(1-fy)*  fz;
        rho(ixloc  , (iy+1)%N, iz      ) += (1-fx)*  fy  *(1-fz);
        rho(ixloc  , (iy+1)%N, (iz+1)%N) += (1-fx)*  fy  *  fz  ;
        rho(ixloc+1, iy      , iz      ) +=   fx  *(1-fy)*(1-fz);
        rho(ixloc+1, iy      , (iz+1)%N) +=   fx  *(1-fy)*  fz  ;
        rho(ixloc+1, (iy+1)%N, iz      ) +=   fx  *  fy  *(1-fz);
        rho(ixloc+1, (iy+1)%N, (iz+1)%N) +=   fx  *  fy  *  fz  ;
    }

    /* Transfer boundary layers across slabs (via a shift operation) */
    int source = (nprocs+me-1) % nprocs;        // process 0 sends to process 1...
    int dest = (me+1) % nprocs;                 // and receives from process nprocs-1
    int nboundary = N*2*(N/2+1);                // number of fftw_real's in boundary layer
    fftw_real* boundary = (fftw_real*)calloc(nboundary, sizeof(fftw_real));  // receive buffer
    MPI_Sendrecv(&rho(nxloc,0,0), nboundary, MPI_DOUBLE, dest, 0,
                 boundary, nboundary, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    for(iy = 0; iy < N; iy++)
        for(iz = 0; iz < N; iz++)
            rho(0,iy,iz) += boundary[iy*(2*(N/2+1)) + iz];
    free(boundary);

    /* Get total number of galaxies across all processes */
//    printf("Process %d: ngals = %d\n", me, ngals);
    MPI_Allreduce(MPI_IN_PLACE, &ngals, 1, MPI_INT, MPI_SUM, comm);
///    printf("After Allreduce: process %d, ngals = %d\n", me, ngals);

    /* Sanity check: grid values should sum to ngals */
    {
        double sum = 0;
        for(ixloc = 0; ixloc < nxloc; ixloc++)
            for(iy = 0; iy < N; iy++)
                for(iz = 0; iz < N; iz++)
                    sum += rho(ixloc,iy,iz);
//        printf("Process %d: sum = %g\n", me, sum);
        if(me == 0) {
            MPI_Reduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            printf("Number density: grid sum = %g, ngals = %d\n", sum, ngals);
        }
        else {
            MPI_Reduce(&sum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        }
    }

    /* Convert to density contrast $\delta(\vec{x})$ */
    double norm = ngals/pow(N, 3);
    for(ixloc = 0; ixloc < nxloc; ixloc++)
        for(iy = 0; iy < N; iy++)
            for(iz = 0; iz < N; iz++)
                rho(ixloc,iy,iz) = rho(ixloc,iy,iz)/norm - 1;

    /* Sanity check: grid values should sum to zero */
    {
        double sum = 0;
        for(ixloc = 0; ixloc < nxloc; ixloc++)
            for(iy = 0; iy < N; iy++)
                for(iz = 0; iz < N; iz++)
                    sum += rho(ixloc,iy,iz);
        if(me == 0) {
            MPI_Reduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            printf("Density contrast: grid sum = %g\n", sum);
        }
        else {
            MPI_Reduce(&sum, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        }
    }

    /* Fourier transform to get $\delta_{\vec{k}}$ */
    rho.fft(false);

    /* Bin amplitudes $|\delta_{\vec{k}}|^2$, correcting for CIC windowing */
    int jxloc, jx, jy, jz, bin;
    double kx, ky, kz, k, pk, mult;
    double kfun = 2*M_PI/L;
    for(jxloc = 0; jxloc < nxloc; jxloc++) {
        jx = ixmin + jxloc;
        kx = kfun * (jx - N*(jx > N/2));
        for(jy = 0; jy < N; jy++) {
            ky = kfun * (jy - N*(jy > N/2));
            for(jz = 0; jz <= N/2; jz++) {
                kz = kfun * jz;
                k = sqrt(kx*kx + ky*ky + kz*kz);
                bin = int(nbins*(k - kmin)/(kmax - kmin));
                if(bin >= 0 && bin < nbins) {
                    mult = 2 - (jz == 0) - (jz == N/2); // multiplicity: 2 for iz != and iz != N/2, since |\delta_{\vec{k}}|^2 = |\delta_{-\vec{k}}|^2
                    pk = pow(rho.c(jxloc,jy,jz).re, 2) + pow(rho.c(jxloc,jy,jz).im, 2);
                    pk *= mult / pow(W(kx*L/N, ky*L/N, kz*L/N), 2);
                    power[bin] += pk;
                    count[bin] += mult;
                }
            }
        }
    }

    /* Accumalate counts across all processes */
    if(me == 0) {
        MPI_Reduce(MPI_IN_PLACE, &power[0], nbins, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(MPI_IN_PLACE, &count[0], nbins, MPI_DOUBLE, MPI_SUM, 0, comm);
    }
    else {
        MPI_Reduce(&power[0], NULL, nbins, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&count[0], NULL, nbins, MPI_DOUBLE, MPI_SUM, 0, comm);
    }

    /* Rescale and output */
    if(me == 0) {
        fprintf(fout, "# Power spectrum computed by gmock\n");
        fprintf(fout, "# N = %d, L = %g, ngals = %d\n", N, L, ngals);
        for(bin = 0; bin < nbins; bin++) {
            k = kmin + (bin + 0.5)*(kmax - kmin)/nbins;
            if(count[bin] > 0)
                pk = pow(L/(N*N), 3) * power[bin]/count[bin] - pow(L, 3)/ngals;
            else
                pk = 0;
            fprintf(fout, "%e %e\n", k, pk);
        }
    }
}
#endif // HAVE_MPI

int main(int argc, char* argv[]) {
    int nprocs = 1, me = 0;     // MPI parameters (default is compatible with non-MPI build)

#ifdef HAVE_MPI
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);
#endif

    /* Parse command line and configuration options */
    Config cfg = cfg_new();
    int opt;
    while((opt = getopt(argc, argv, "hc:")) != -1) {
        switch(opt) {
        case 'h':
            printf("Usage: gmock -c cfgfile\n");
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fprintf(stderr, "Invalid command line options.  Try 'gmock -h'.\n");
            return 1;
        }
    }
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Debugging */
    if(me == 0) {
        printf("Config options:\n");
        cfg_write(cfg, stdout);
        printf("\n");
    }

    /* Quick fix: grab all config entries starting with "gmock." */
    {
        Config gmockcfg = cfg_new_sub(cfg, "gmock.", 1);
        cfg_destroy(cfg);
        cfg = gmockcfg;
    }

//    if(!cfg_has_keys(cfg, "pkfile,N,L,f,b,rhobar,sigma,outfile", ",")) {
    if(!cfg_has_keys(cfg, "pkfile,N,L,observer,f,b,rhobar,Rsmooth,outfile", ",")) {
        fprintf(stderr, "Missing configuration options.\n");
        return 1;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    int N = cfg_get_int(cfg, "N");
    double L = cfg_get_double(cfg, "L");
    double observer[4] = { 0, 0, 0, 1 };       // default to observer at the origin
    cfg_get_array_double(cfg, "observer", 4, observer);
    double f = cfg_get_double(cfg, "f");
    double b = cfg_get_double(cfg, "b");
    double rhobar = cfg_get_double(cfg, "rhobar");
//    double sigma = cfg_get_double(cfg, "sigma");
    double Rsmooth = cfg_get_double(cfg, "Rsmooth");
    const char* outfile = cfg_get(cfg, "outfile");

    /* Impose parameter restrictions: things are way easier this way */
    if((N % nprocs) != 0) {
        fprintf(stderr, "Error: N must be a multiple of nprocs (N = %d, nprocs = %d)\n", N, nprocs);
        return 1;
    }

    /* Initialize slab information */
    Slab slab;
    slab.N = N;
    slab.L = L;
    slab.nxloc = N/nprocs;
    slab.ny = N;
    slab.nz = N;
    slab.ixmin = me*slab.nxloc;
    Grid& rho = slab.rho;
    Grid& vx = slab.vx;
    Grid& vy = slab.vy;
    Grid& vz = slab.vz;
    int nxloc = slab.nxloc;

    /* Make sure we can open output file */
    FILE* fout = NULL;
    if(me == 0) {
        fout = fopen(outfile, "w");
        if(!fout) {
            fprintf(stderr, "gmock: cannot open '%s' for writing.\n", outfile);
            return 1;
        }
    }

    /* Seed random number generator.  All processes choose the same seed: this
     * is important for generating the density field. */
    unsigned long seed = cfg_has_key(cfg, "seed") ? cfg_get_uint(cfg, "seed") : DEFAULT_SEED;
    rng_init(seed);

    if(me == 0)
        printf("Loading power spectrum\n");

    /* Load power spectrum from file */
    Spline P = LinearSpline(pkfile);

    if(me == 0)
        printf("Preparing grid\n");

    /* Allocate grids */
#ifdef HAVE_MPI
    rho.initialize(comm, N, N, N);
    vx.initialize(comm, N, N, N);
    vy.initialize(comm, N, N, N);
    vz.initialize(comm, N, N, N);
#else
    rho.initialize(N, N, N);
    vx.initialize(N, N, N);
    vy.initialize(N, N, N);
    vz.initialize(N, N, N);
#endif // HAVE_MPI

    if(me == 0)
        printf("Generating density and velocity fields\n");

    generate_fields(N, L, f, P, Rsmooth, rho, vx, vy, vz);

#if 0
    /* Rescale density contrast to have desired sigma */
    {
        double var = 0.;
        for(int ixloc = 0; ixloc < nxloc; ixloc++)
            for(int iy = 0; iy < N; iy++)
                for(int iz = 0; iz < N; iz++)
                    var += pow(rho(ixloc,iy,iz), 2);

#ifdef HAVE_MPI
        if(me == 0)
            MPI_Reduce(MPI_IN_PLACE, &var, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        else
            MPI_Reduce(&var, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

        double norm = sigma*sqrt(N*N*N/var);
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, comm);
        if(me == 0)
            printf("  norm = %g\n", norm);
#else
        double norm = sigma*sqrt(N*N*N/var);
#endif

        for(int ixloc = 0; ixloc < nxloc; ixloc++) {
            for(int iy = 0; iy < N; iy++) {
                for(int iz = 0; iz < N; iz++) {
                    rho(ixloc,iy,iz) *= norm;
                    vx(ixloc,iy,iz) *= norm;
                    vy(ixloc,iy,iz) *= norm;
                    vz(ixloc,iy,iz) *= norm;
                }
            }
        }
    }
#endif

    /* Convert density contrast to proper number density */
    double numneg = 0;
    for(int ixloc = 0; ixloc < nxloc; ixloc++) {
        for(int iy = 0; iy < N; iy++) {
            for(int iz = 0; iz < N; iz++) {
                if(rho(ixloc,iy,iz) < -1/b)
                    numneg += 1;
                rho(ixloc,iy,iz) = rhobar * fmax(0, 1 + b*rho(ixloc,iy,iz));
            }
        }
    }
    printf("Process %d: density field %g percent negative\n", me, 100*numneg/(nxloc*N*N));

    if(me == 0)
        printf("Transferring boundary layers\n");

    /* Use the extra memory allocated by Grid (nxloc+1 instead of nxloc) to
     * store the grid values from the adjacent process, process me+1.  At the
     * same time, we need to send our boundary grid values to process me-1. */
    int dest = (nprocs+me-1) % nprocs;
    int source = (me+1) % nprocs;
    int count = N*2*(N/2+1);
    assert(sizeof(fftw_real) == sizeof(double));
    MPI_Sendrecv(&rho(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &rho(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vx(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vx(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vy(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vy(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vz(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vz(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);

    if(me == 0)
        printf("Sampling galaxies from density field\n");

    /* Poisson sample from density field */
    vector<Galaxy> local_galaxies;
    int num = sample_galaxies(slab, observer, local_galaxies);

#ifdef HAVE_MPI
    /* For fun, calculate the power spectrum for the galaxies just sampled. */
    ComputePowerSpectrum(comm, local_galaxies, slab, "pk-gmock.dat");
#endif

    if(me == 0)
        printf("Writing galaxies to '%s'\n", outfile);

    /* Gather galaxies and write to file */
#ifdef HAVE_MPI
    /* Determine the total number of galaxies, as well as the maximum number
     * of galaxies owned by any single process */
    int totalnum = 0, maxnum = 0;
    MPI_Reduce(&num, &totalnum, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&num, &maxnum, 1, MPI_INT, MPI_MAX, 0, comm);

    /* Prepare galaxy output file header */
    if(me == 0) {
        Config opts = cfg_new();
        cfg_set(opts, "coordsys", "cartesian");
        abn_write_header(fout, totalnum, "4f", opts);
    }

    /* Write root process's galaxies to file */
    if(me == 0) {
        fwrite(&local_galaxies[0], sizeof(float), 4*local_galaxies.size(), fout);
        local_galaxies.clear(); // free up memory
    }

    /* Allocate memory so root process can hold maxnum galaxies */
    vector<Galaxy> pbuf;
    if(me == 0)
        pbuf.resize(maxnum);

    /* Send galaxies to root process and write to file, one process at a time */
    for(int i = 1; i < nprocs; i++) {
        if(i == me) {
//            printf("Process %d sending %d galaxies to root\n", i, num);
            MPI_Send(&local_galaxies[0], 4*num, MPI_FLOAT, 0, 0, comm);
        }
        if(me == 0) {
            int count;
            MPI_Status status;
            MPI_Recv(&pbuf[0], 4*maxnum, MPI_FLOAT, i, 0, comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &count);
//            printf("Received %d galaxies from process %d\n", count/3, i);
            fwrite(&pbuf[0], sizeof(float), count, fout);
        }
        MPI_Barrier(comm);
    }
#else
    Config opts = cfg_new();
    cfg_set(opts, "coordsys", "cartesian");
    abn_write(fout, &local_galaxies[0], num, "4f", opts);
#endif // HAVE_MPI

    /* Clean up */
    if(me == 0)
        fclose(fout);
    rho.cleanup();
    vx.cleanup();
    vy.cleanup();
    vz.cleanup();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}
