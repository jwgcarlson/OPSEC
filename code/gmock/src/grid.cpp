#include <cassert>
#include <cstring>

#include "grid.h"

Grid::Grid() {
    nx = ny = nz = nxloc = ixmin = nyloc_t = iymin_t = local_size = 0;
    data = NULL;
    fft_plan = ifft_plan = NULL;
}

Grid::~Grid() {
    cleanup();
}

#ifdef HAVE_MPI
void Grid::initialize(MPI_Comm comm, int nx_, int ny_, int nz_) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    fft_plan = rfftw3d_mpi_create_plan(comm, nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    ifft_plan = rfftw3d_mpi_create_plan(comm, nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    rfftwnd_mpi_local_sizes(fft_plan, &nxloc, &ixmin, &nyloc_t, &iymin_t, &local_size);
#else
void Grid::initialize(int nx_, int ny_, int nz_) {
    nx = nx_;
    ny = ny_;
    nz = nz_;
    rfftwnd_plan plan = rfftw3d_create_plan(nx, ny, nz, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    rfftwnd_plan iplan = rfftw3d_create_plan(nx, ny, nz, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
    nxloc = nx;
    nyloc_t = ny;
    ixmin = iymin_t = 0;
    local_size = nx * ny * 2*(nz/2+1);
#endif // HAVE_MPI

    /* Allocate extra storage so that each process can hold the boundary
     * layer from the adjacent process */
    assert(local_size == nxloc*ny*2*(nz/2+1));
    local_size = (nxloc+1)*ny*2*(nz/2+1);
    data = (fftw_real*) malloc(local_size*sizeof(fftw_real));
}

void Grid::cleanup() {
#ifdef HAVE_MPI
    rfftwnd_mpi_destroy_plan(fft_plan);
    rfftwnd_mpi_destroy_plan(ifft_plan);
#else
    rfftwnd_destroy_plan(fft_plan);
    rfftwnd_destroy_plan(ifft_plan);
#endif // HAVE_MPI
    free(data);
    nx = ny = nz = nxloc = ixmin = nyloc_t = iymin_t = local_size = 0;
    data = NULL;
    fft_plan = ifft_plan = NULL;
}

#ifdef HAVE_MPI
void Grid::fft(bool transpose) {
    fftwnd_mpi_output_order order = transpose ? FFTW_TRANSPOSED_ORDER : FFTW_NORMAL_ORDER;
    rfftwnd_mpi(fft_plan, 1, data, NULL, order);
}

void Grid::zero() {
    memset(data, 0, local_size*sizeof(fftw_real));
}

void Grid::ifft(bool transpose) {
    fftwnd_mpi_output_order order = transpose ? FFTW_TRANSPOSED_ORDER : FFTW_NORMAL_ORDER;
    rfftwnd_mpi(ifft_plan, 1, data, NULL, order);
}
#else
void Grid::fft(bool) {
    rfftwnd_one_real_to_complex(fft_plan, data, NULL);
}

void Grid::ifft(bool) {
    rfftwnd_one_complex_to_real(ifft_plan, data, NULL);
}
#endif // HAVE_MPI
