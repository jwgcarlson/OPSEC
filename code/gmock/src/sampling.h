#ifndef SAMPLING_H
#define SAMPLING_h

#include <vector>

#include "grid.h"
#include "particle.h"
#include "Spline.h"

/* Sample x \in [0,1] from the linear 1-D distribution
 *   f(x) = a0*(1-x) + a1*x */
double sample1d(double a0, double a1);

/* Sample (x,y,z) \in [0,1]^3 from a 3-D distribution f(x,y,z) that is
 * tri-linear in (x,y,z), with the values at the corners given by a_{ijk} */
void sample3d(double a000, double a001, double a010, double a011,
              double a100, double a101, double a110, double a111,
              double& x, double& y, double& z);

/* Generate random Fourier modes $\delta_{\vec{k}}$ according to the given
 * P(k), then Fourier transform to obtain the real space fields.  The velocity
 * field (vx,vy,vz) is given in Fourier space by
 *   \vec{v}_{\vec{k}} = (i\vec{k}/k^2) f \delta_{\vec{k}}
 * The density and velocity fields are smoothed on the scale R. */
void generate_fields(int N, double L, double f, Spline P, double R, Grid& delta, Grid& vx, Grid& vy, Grid& vz);

#if 0
/* Inverse CIC: return the value of the field at the point (x,y,z) \in [0,1]^3 */
double icic(double a000, double a001, double a010, double a011,
            double a100, double a101, double a110, double a111,
            double x, double y, double z);

/* Sample particles.
 *   x = x0 + (ix+fx)*dL
 *   y = y0 + (iy+fy)*dL
 *   z = z0 + (iz+fz)*dL
 * Return the number of particles added. */
int sample_particles(Grid& rho,                         /* number density field */
                     Grid& vx, Grid& vy, Grid& vz,      /* velocity field */
                     double dL,                         /* side length of one grid cell */
                     int ixbegin, int ixend, int nx,
                     int iybegin, int iyend, int ny,
                     int izbegin, int izend, int nz,
                     double x0, double y0, double z0,   /* offset for particles */
                     std::vector<Particle>& particles   /* array of particles to add to */
                    );
#endif

#endif // SAMPLING_H
