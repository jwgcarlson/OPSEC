#ifndef GRID_H
#define GRID_H

#include <cassert>
#include <rfftw.h>

#ifdef OPSEC_USE_MPI
/* MPI Grid, distributed as slabs among processes */

#include <rfftw_mpi.h>

struct Grid {
    int nx, ny, nz;             // dimensions of full grid
    int nxloc, ixmin;           // size and starting index of local slab
    int nyloc_t, iymin_t;       // size and starting index of transposed local slab
    int local_size;             // total memory allocated for local data
    fftw_real* data;            // raw local data
    rfftwnd_mpi_plan fft_plan;  // FFTW2 plan for forward FFT
    rfftwnd_mpi_plan ifft_plan; // FFTW2 plan for inverse FFT

    Grid();
    ~Grid();
    void initialize(MPI_Comm comm, int nx, int ny, int nz);
    void cleanup();
    void zero();
    void fft(bool transpose = false);
    void ifft(bool transpose = false);

    fftw_real& r(int ixloc, int iy, int iz) const {
        return data[(ixloc*ny + iy)*(2*(nz/2+1)) + iz];
    }

    fftw_complex& c(int ixloc, int iy, int iz) const {
        return ((fftw_complex*)data)[(ixloc*ny + iy)*(nz/2+1) + iz];
    }

    fftw_real& r_t(int ix, int iyloc, int iz) const {
        return data[(iyloc*nx + ix)*(2*(nz/2+1)) + iz];
    }

    fftw_complex& c_t(int ix, int iyloc, int iz) const {
        return ((fftw_complex*)data)[(iyloc*nx + ix)*(nz/2+1) + iz];
    }

    /* Shorthand for most common case */
    fftw_real& operator()(int ixloc, int iy, int iz) const {
        return r(ixloc,iy,iz);
    }
};

#else
/* Non-MPI grid, owned by a single process */

struct Grid {
    int nx, ny, nz;             // dimensions of full grid
    int nxloc, ixmin;           // size and starting index of local slab
    int nyloc_t, iymin_t;       // size and starting index of transposed local slab
    int local_size;             // total memory allocated for local data
    fftw_real* data;            // raw local data
    rfftwnd_plan fft_plan;      // FFTW2 plan for forward FFT
    rfftwnd_plan ifft_plan;     // FFTW2 plan for inverse FFT

    Grid();
    ~Grid();
    void initialize(int nx, int ny, int nz);
    void cleanup();
    void zero();
    void fft(bool);
    void ifft(bool);

    fftw_real& r(int ixloc, int iy, int iz) const {
        return data[(ixloc*ny + iy)*(2*(nz/2+1)) + iz];
    }

    fftw_complex& c(int ixloc, int iy, int iz) const {
        return ((fftw_complex*)data)[(ixloc*ny + iy)*(nz/2+1) + iz];
    }

    /* For the non-MPI case, we never actually transpose the grid. */
    fftw_real& r_t(int ix, int iyloc, int iz) const {
        return data[(ix*ny + iyloc)*(2*(nz/2+1)) + iz];
    }

    fftw_complex& c_t(int ix, int iyloc, int iz) const {
        return ((fftw_complex*)data)[(ix*ny + iyloc)*(nz/2+1) + iz];
    }

    /* Shorthand for most common case */
    fftw_real& operator()(int ixloc, int iy, int iz) const {
        return r(ixloc,iy,iz);
    }
};

#endif // OPSEC_USE_MPI

#endif // GRID_H
