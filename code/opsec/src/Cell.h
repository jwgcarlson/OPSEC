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


/* Coordinate system options. */
enum {
    CoordSysUndefined = 0,
    CoordSysCartesian = 1,
    CoordSysSpherical = 2
};

/* A point in either Cartesian or spherical coordinates. */
struct Point {
    union { double x1; double x; double r; };
    union { double x2; double y; double mu; };
    union { double x3; double z; double phi; };
};


/* Structure representing a non-empty cell within the survey grid. */
struct Cell {
    /* Cell number within list of non-empty cells: 0 <= a < Ncells */
    int a;

    /* Cell index within survey grid: 0 <= d1 < N1, 0 <= d2 < N2, 0 <= d3 < N3 */
    union { int d1; int dx; int dr; };
    union { int d2; int dy; int dmu; };
    union { int d3; int dz; int dphi; };

    /* Cell bounds */
    union { double min1; double xmin; double rmin; };
    union { double max1; double xmax; double rmax; };
    union { double min2; double ymin; double mumin; };
    union { double max2; double ymax; double mumax; };
    union { double min3; double zmin; double phimin; };
    union { double max3; double zmax; double phimax; };

    /* Effective cell volume */
    double Veff;

    /* Expected number of galaxies within cell: $\bar{N} = \int_V \bar{n} dV$ */
    double Nbar;
};

/* ABN format string for cells. */
#define CELL_FMT_STRING "4i8d"

Cell* ReadCells(const char* cellfile, int* Ncells);

#endif // CELL_H
