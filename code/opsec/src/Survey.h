#ifndef SURVEY_H
#define SURVEY_H

#include <cstdio>
#include <vector>

#include "Cell.h"
#include "cfg.h"

class SelectionFunc;
class SeparationFunc;

/* Data structure representing a galaxy in either Cartesian or spherical coordinates. */
struct Galaxy {
    union { float x1; float x; float r; };
    union { float x2; float y; float mu; };
    union { float x3; float z; float phi; };
    float w;    // galaxy weight
};

/* Survey
 * 
 * Abstract represention of a galaxy survey, including a separation function
 * (function for computing the separation between pairs of galaxies), the
 * selection function (expected number density of galaxies) and the observed
 * galaxy positions (as a list of Galaxy objects). */
class Survey {
public:
    Survey(int coordsys, Config cfg = NULL);
    virtual~ Survey();

    /* Return the coordinate system for this survey. */
    int GetCoordinateSystem() const;

    /* Return an empty cell object for the given grid position. */
    Cell CreateEmptyCell(int d1, int d2, int d3) const;

    /* Return the grid index G = (d*N2 + e)*N3 + f of the cell containing the
     * given point, or -1 if the point does not lie within the survey region. */
    int GetGridIndex(double x1, double x2, double x3) const;

    /* Return a SeparationFunc object for this survey. */
    virtual SeparationFunc GetSeparationFunction() = 0;

    /* Return a callable SelectionFunc object representing the selection
     * function (expected number density of galaxies) for this survey. */
    virtual SelectionFunc GetSelectionFunction() = 0;

    /* Populate the vector 'gals' with the galaxies in the survey.  If the
     * supplied vector is non-empty, the galaxies are appended. */
    virtual void GetGalaxies(std::vector<Galaxy>& gals) = 0;

    /* Return a configuration dictionary with all relevant options. */
    virtual Config GetConfigurationOptions() const;

public:
    /* Coordinate system flag, either CoordSysCartesian or CoordSysSpherical */
    int coordsys;

    /* Coordinate-aligned survey bounding region */
    union { double Min1; double XMin; double RMin; };
    union { double Max1; double XMax; double RMax; };
    union { double Min2; double YMin; double MuMin; };
    union { double Max2; double YMax; double MuMax; };
    union { double Min3; double ZMin; double PhiMin; };
    union { double Max3; double ZMax; double PhiMax; };

    /* Bounding region subdivision */
    union { int N1; int Nx; int Nr; };
    union { int N2; int Ny; int Nmu; };
    union { int N3; int Nz; int Nphi; };
};

/* Load a concrete Survey instance based on the configuration options
 * specified in 'cfg'.  The configuration option "survey" must be set, and
 * must correspond to a fully implemented subclass of Survey. */
Survey* InitializeSurvey(Config cfg);

#endif // SURVEY_H
