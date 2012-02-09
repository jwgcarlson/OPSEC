#ifndef SIG_H
#define SIG_H

#include <cstdio>

#include "Cell.h"
#include "Model.h"
#include "opsec.h"

/* Compute the local nloc-by-n block of the n-by-n signal matrix, starting at
 * row amin.  Assume spherical coordinates, and use azimuthal symmetry to
 * minimize computation. */
real* ComputeSignalMatrixS(size_t n, size_t nloc, int amin, const Cell* cells,
                           int Nr, int Nmu, int Nphi, const XiFunc& xi,
                           double epsrel = 1e-5, double epsabs = 1e-10);

/* Compute the local nloc-by-n block of the n-by-n signal matrix, starting at
 * row amin.  Assume Cartesian coordinates and translational symmetry (real
 * space, not redshift space) to minimize computation. */
real* ComputeSignalMatrixC(size_t n, size_t nloc, int amin, const Cell* cells,
                           int Nx, int Ny, int Nz, const XiFunc& xi,
                           double epsrel = 1e-5, double epsabs = 1e-10);

/* Read the header from the modes.abn file.  Return the file pointer to the
 * opened file, set to point at the first element of the first mode vector.
 * If non-NULL, set Nmodes and/or Ncells. */
FILE* ReadModesHeader(const char* modefile, int* Nmodes = NULL, int* Ncells = NULL);

/* Read KL mode coefficients from file. */
real* ReadModes(const char* modefile, int* Nmodes = NULL, int* Ncells = NULL);

/* Auxiliary routines, should not be public anymore.
 * Except yes they have to be for dot to be able to compute variances in counts. */
double ComputeSignalS(const Cell& c1, const Cell& c2, const XiFunc& xi, double epsrel = 1e-5, double epsabs = 1e-10);
double ComputeSignalC(const Cell& c1, const Cell& c2, const XiFunc& xi, double epsrel = 1e-5, double epsabs = 1e-10);

#endif // SIG_H
