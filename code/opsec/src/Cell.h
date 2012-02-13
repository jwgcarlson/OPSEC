#ifndef CELL_H
#define CELL_H

/* The survey volume is surrounded by a regular region B defined by spherical
 * coordinate boundaries:
 *   B = { (r,theta,phi) | RMin < r < RMax, MuMin < mu < MuMax, PhiMin < phi < PhiMax }
 * This bounding region is then partitioned along the r, mu, and phi directions
 * into a coordinate-aligned grid, with
 *   Nr radial divisions
 *   Nmu polar angle divisions
 *   Nphi azimuthal angle divisions
 *
 * Each cell within this partitioning is identified uniquely by integers
 * (d,e,f) with
 *   0 <= d < Nr, 0 <= e < Nmu, 0 <= f < Nphi,
 * and the boundaries of each individual cell are
 *   rmin = RMin +     d*(RMax - RMin)/Nr
 *   rmax = RMin + (d+1)*(RMax - RMin)/Nr
 *   mumin = MuMin +     e*(MuMax - MuMin)/NMu
 *   mumax = MuMin + (e+1)*(MuMax - MuMin)/NMu
 *   phimin = ThetaMin +     f*(PhiMax - PhiMin)/Nphi
 *   phimax = ThetaMin + (f+1)*(PhiMax - PhiMin)/Nphi
 * The cell volume may be computed easily in terms of the coordinate bounds:
 *   V = (rmax^3 - rmin^3)/3 * (mumax - mumin) * (phimax - phimin).
 * The triplet (d,e,f) is compressed into the single integer 'G', defined by
 *   G = f + Nphi*e + Nmu*Nphi*d;
 * This can be easily inverted (exploiting C rules for integer arithmetic) to
 * obtain (d,e,f) as
 *   d = G / (Nmu*Nphi);
 *   e = (G / Nphi) % Nmu;
 *   f = G % Nphi
 *
 * Some cells within the grid B may not be contained in the survey volume.
 * These cells are discarded, so the number of cells may be less than
 * Nr*Nmu*Nphi.  The cell number 'a' is a consecutive index that uniquely
 * labels each cell.
 *
 * The survey volume may intersect only a part of a given cell, with the
 * selection function vanishing over some fraction fzero of the cell.  That
 * cell is given an effective volume Veff = (1-fzero)*V, with V computed from
 * the coordinate bounds.  The expected number of galaxies in the cell (in
 * the absence of clustering) is defined by the integral
 *   Nbar = \int_V \bar{n}(\vec{s}) d^3s
 * over the whole volume V of the cell.
 *
 * Corresponding to each cell is a basis function defined to have the constant
 * value 1/sqrt(Nbar) within the cell, and zero elsewhere.  This value is
 * chosen to pre-whiten the noise matrix (i.e. make it equal the identity).
 *
 * The above discussion refers to spherically regular cells, i.e. rectangular
 * volumes in (r,mu,phi) space.  For Cartesian cells, replace (r,mu,phi) with
 * (x,y,z). */

/* NB: Be sure to update MPI_Datatype definition in ReadCells() if the Cell
 * structure changes! */


/* Coordinate system options */
enum {
    CoordSysCartesian = 0,
    CoordSysSpherical = 1
};


/* Abuse anonymous unions to allow the same struct to represent either a
 * spherical cell or a Cartesian cell. */
struct Cell {
    int a;                      // cell number (0 <= a < Ncells)
    int G;                      // cell location within grid (0 <= G < N1*N2*N3)
    union { double min1; double xmin; double rmin; };
    union { double max1; double xmax; double rmax; };
    union { double min2; double ymin; double mumin; };
    union { double max2; double ymax; double mumax; };
    union { double min3; double zmin; double phimin; };
    union { double max3; double zmax; double phimax; };
    double Veff;                // effective volume of cell
    double Nbar;                // expected number of galaxies within cell (Nbar = \int_V \bar{n} dV)
};

Cell* ReadCells(const char* cellfile, int* Ncells);

#endif // CELL_H
