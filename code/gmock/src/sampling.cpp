#include <cassert>
#include <cmath>
#include <cstdio>

#include <vector>
using std::vector;

#include "rng.h"
#include "sampling.h"

double sample1d(double a0, double a1) {
    assert(a0 >= 0 && a1 >= 0);
    double norm = 0.5*(a0 + a1);
    a0 /= norm;  a1 /= norm;

    double u = rng_uniform();
    if(a0 == a1)
        return u;
    else
        return (-a0 + sqrt(a0*a0 + 2*(a1-a0)*u))/(a1 - a0);
}

void sample3d(double a000, double a001, double a010, double a011,
              double a100, double a101, double a110, double a111,
              double& x, double& y, double& z)
{
    assert(a000 >= 0 && a001 >= 0 && a010 >= 0 && a011 >= 0 && a100 >= 0 && a101 >= 0 && a110 >= 0 && a111 >= 0);

    /* Normalize coefficients */
    double norm = 0.125*(a000 + a001 + a010 + a011 + a100 + a101 + a110 + a111);
    a000 /= norm; a001 /= norm; a010 /= norm; a011 /= norm; a100 /= norm; a101 /= norm; a110 /= norm; a111 /= norm;

    double a00 = 0.5*(a000 + a001);
    double a01 = 0.5*(a010 + a011);
    double a10 = 0.5*(a100 + a101);
    double a11 = 0.5*(a110 + a111);

    /* Choose x position from 1-D marginal distribution, f_X(x) = a0*(1-x) + a1*x */
    double a0 = 0.5*(a00 + a01);
    double a1 = 0.5*(a10 + a11);
    x = sample1d(a0, a1);

    /* Choose y position from conditional distribution,  f_{Y|X}(y|x) = b0*(1-y) + b1*y */
    double b0 = a00*(1-x) + a10*x;
    double b1 = a01*(1-x) + a11*x;
    y = sample1d(b0, b1);

    /* Choose z position form conditional distribution, f_{Z|X,Y}(z|x,y) = c0*(1-z) + c1*z */
    double c0 = a000*(1-x)*(1-y) + a010*(1-x)*y + a100*x*(1-y) + a110*x*y;
    double c1 = a001*(1-x)*(1-y) + a011*(1-x)*y + a101*x*(1-y) + a111*x*y;
    z = sample1d(c0, c1);
}

void generate_fields(int N, double L, double f, Spline P, double R, Grid& delta, Grid& vx, Grid& vy, Grid& vz) {
    double V = L*L*L;
    double damping, R2 = R*R;
    int iymin_t = delta.iymin_t, nyloc_t = delta.nyloc_t;
    int ix, iy, iz;             // mode index
    int cix, ciy;               // conjugate mode index
    int iyloc, ciyloc;          // mode index within current slab
    double kx, ky, kz, k, k2;
    double kfun = 2*M_PI/L;     // fundamental frequency
    double knyq = M_PI*N/L;     // Nyquist frequency
    double u1, u2;              // uniform random variates
    double mod, arg, re, im;
    for(ix = 0; ix < N; ix++) {
        cix = (N - ix) % N;
        kx = kfun * (ix - N*(ix > N/2));
        for(iy = 0; iy < N; iy++) {
            ciy = (N - iy) % N;
            iyloc = iy - iymin_t;
            ciyloc = ciy - iymin_t;
            ky = kfun * (iy - N*(iy > N/2));
            for(iz = 0; iz <= N/2; iz++) {
                kz = kfun * iz;

                k2 = kx*kx + ky*ky + kz*kz;
                k = sqrt(k2);

                /* Sample the random number generator, even if the mode doesn't
                 * belong to our slab.  This is done to ensure that reality
                 * conditions (delta_k* = delta_{-k}) are satisfied across
                 * processes. */
                u1 = rng_uniform();
                u2 = rng_uniform();

                if(k > knyq) {
                    /* Cut off power isotropically above the Nyquist frequency */
                    re = im = 0.;
                }
                else {
                    mod = sqrt(-P(k)/V * log(1-u1));
                    arg = 2*M_PI * u2;
                    re = mod*cos(arg);      // Re \delta_{\vec{k}}
                    im = mod*sin(arg);      // Im \delta_{\vec{k}}
                }

                /* With FFTW, the complex Fourier coefficients for a real-valued
                 * density field are arranged in a grid of dimension
                 * N*N*(N/2+1).  This grid contains 2N^2 redundant real values,
                 * which are determined from the other coefficients by the
                 * reality condition
                 *   \delta_{-\vec{k}} = \delta_{\vec{k}}^*
                 * The redundant coefficients can be characterized by
                 *   ((iz == 0 || iz == N/2) && iy > N/2) ||
                 *   ((iz == 0 || iz == N/2) && (iy == 0 || iy == N/2) && ix > N/2)
                 * This amounts to
                 *   2*(N/2-1)*N + 2*2*(N/2-1) = N^2 - 4
                 * completely redundant complex coefficients.  Moreover the 8
                 * coefficients with
                 *   (iz == 0 || iz == N/2) && (iy == 0 || iy == N/2) && (ix == 0 || ix == N/2)
                 * must be real, which combined with the above accounts for the
                 * 2N^2 redundant real values. */

                if((iz == 0 || iz == N/2) && (iy > N/2 || ((iy == 0 || iy == N/2) && ix > N/2)))
                    continue;   // skip redundant mode

                if((iz == 0 || iz == N/2) && (iy == 0 || iy == N/2) && (ix == 0 || ix == N/2))
                    im = 0.;    // real-valued coefficient

                if(k2 == 0) {
                    re = im = 0;        // zero mode
                    k2 = 1.;            // prevent nan's in velocity coefficients
                }

                /* Apply smoothing to the density and velocity fields */
                damping = exp(-k2*R2/2);
                re *= damping;
                im *= damping;

                /* If the mode belongs to our slab, set the appropriate local
                 * grid elements */
                if(iyloc >= 0 && iyloc < nyloc_t) {
                    delta.c_t(ix,iyloc,iz).re = re;
                    delta.c_t(ix,iyloc,iz).im = im;
                    vx.c_t(ix,iyloc,iz).re = -kx/k2 * f * im;
                    vx.c_t(ix,iyloc,iz).im = kx/k2 * f * re;
                    vy.c_t(ix,iyloc,iz).re = -ky/k2 * f * im;
                    vy.c_t(ix,iyloc,iz).im = ky/k2 * f * re;
                    vz.c_t(ix,iyloc,iz).re = -kz/k2 * f * im;
                    vz.c_t(ix,iyloc,iz).im = kz/k2 * f * re;
                }

                /* If the conjugate mode belongs to our slab, set the
                 * appropriate grid elements */
                if((iz == 0 || iz == N/2) && ciyloc >= 0 && ciyloc < nyloc_t) {
                    delta.c_t(cix,ciyloc,iz).re = re;
                    delta.c_t(cix,ciyloc,iz).im = -im;
                    vx.c_t(cix,ciyloc,iz).re = -kx/k2 * f * im;
                    vx.c_t(cix,ciyloc,iz).im = -kx/k2 * f * re;
                    vy.c_t(cix,ciyloc,iz).re = -ky/k2 * f * im;
                    vy.c_t(cix,ciyloc,iz).im = -ky/k2 * f * re;
                    vz.c_t(cix,ciyloc,iz).re = -kz/k2 * f * im;
                    vz.c_t(cix,ciyloc,iz).im = -kz/k2 * f * re;
                }
            }
        }
    }

    /* FFT to obtain real space fields (skip extra transpose step, since we
     * initialized the Fourier coefficients in pre-transposed form) */
    delta.ifft(true);
    vx.ifft(true);
    vy.ifft(true);
    vz.ifft(true);
}


#if 0
double icic(double a000, double a001, double a010, double a011,
            double a100, double a101, double a110, double a111,
            double x, double y, double z)
{
    return a000*(1-x)*(1-y)*(1-z) + a001*(1-x)*(1-y)*z + a010*(1-x)*y*(1-z) + a011*(1-x)*y*z
         + a100*x*(1-y)*(1-z)     + a101*x*(1-y)*z     + a110*x*y*(1-z)     + a111*x*y*z;
}

int sample_particles(Grid& rho,
                     Grid& vx, Grid& vy, Grid& vz,
                     double dL,
                     int ixbegin, int ixend, int nx,
                     int iybegin, int iyend, int ny,
                     int izbegin, int izend, int nz,
                     double x0, double y0, double z0,
                     vector<Particle>& particles)
{
    double dV = dL*dL*dL;       // grid cell volume
    int ntotal = 0;             // number of new particles added

    /* Iterate over grid cells within this chaining mesh cell */
    for(int ix = ixbegin; ix < ixend; ix++) {
        for(int iy = iybegin; iy < iyend; iy++) {
            for(int iz = izbegin; iz < izend; iz++) {
                /* Compute mean density over the grid cell */
                double a000 = rho(ix,        iy,        iz);
                double a001 = rho(ix,        iy,        (iz+1)%nz);
                double a010 = rho(ix,        (iy+1)%ny, iz);
                double a011 = rho(ix,        (iy+1)%ny, (iz+1)%nz);
                double a100 = rho((ix+1)%nx, iy,        iz);
                double a101 = rho((ix+1)%nx, iy,        (iz+1)%nz);
                double a110 = rho((ix+1)%nx, (iy+1)%ny, iz);
                double a111 = rho((ix+1)%nx, (iy+1)%ny, (iz+1)%nz);
                double rhobar = (a000 + a001 + a010 + a011 + a100 + a101 + a110 + a111)/8;

                /* Decide how many new particles to add to this volume element */
                int numnew = (int)rng_poisson(rhobar*dV);
                double fx, fy, fz;      // position of particle within grid cell, (fx,fy,fz) \in [0,1]^3
                Particle p;
                for(int i = 0; i < numnew; i++) {
                    /* Choose where within the grid cell to place the particle */
                    sample3d(a000, a001, a010, a011, a100, a101, a110, a111, fx, fy, fz);
                    p.x = x0 + (ix + fx)*dL;
                    p.y = y0 + (iy + fy)*dL;
                    p.z = z0 + (iz + fz)*dL;

                    /* Compute the velocity at that position */
                    p.vx = icic(vx(ix,iy,iz), vx(ix,iy,(iz+1)%nz), vx(ix,(iy+1)%ny,iz), vx(ix,(iy+1)%ny,(iz+1)%nz), vx((ix+1)%nx,iy,iz), vx((ix+1)%nx,iy,(iz+1)%nz), vx((ix+1)%nx,(iy+1)%ny,iz), vx((ix+1)%nx,(iy+1)%ny,(iz+1)%nz), fx, fy, fz);
                    p.vy = icic(vy(ix,iy,iz), vy(ix,iy,(iz+1)%nz), vy(ix,(iy+1)%ny,iz), vy(ix,(iy+1)%ny,(iz+1)%nz), vy((ix+1)%nx,iy,iz), vy((ix+1)%nx,iy,(iz+1)%nz), vy((ix+1)%nx,(iy+1)%ny,iz), vy((ix+1)%nx,(iy+1)%ny,(iz+1)%nz), fx, fy, fz);
                    p.vz = icic(vz(ix,iy,iz), vz(ix,iy,(iz+1)%nz), vz(ix,(iy+1)%ny,iz), vz(ix,(iy+1)%ny,(iz+1)%nz), vz((ix+1)%nx,iy,iz), vz((ix+1)%nx,iy,(iz+1)%nz), vz((ix+1)%nx,(iy+1)%ny,iz), vz((ix+1)%nx,(iy+1)%ny,(iz+1)%nz), fx, fy, fz);

                    /* Add the particle */
                    particles.push_back(p);
                    ntotal++;
                }
            }
        }
    }

    return ntotal;
}
#endif
