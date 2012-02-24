#ifndef SURVEY_H
#define SURVEY_H

#include <vector>

#include "Cell.h"

class SelectionFunc;
class SeparationFunc;

/* Data structure representing a galaxy in either Cartesian or spherical coordinates. */
struct Galaxy {
    union { float x; float r; };
    union { float y; float mu; };
    union { float z; float phi; };
    float w;    // galaxy weight
};

/* Survey
 * 
 * Abstract represention of a galaxy survey, including a separation function
 * (function for computing the separation between pairs of galaxies), the
 * selection function (expected number density of galaxies) and the observed
 * galaxy positions (as a sequence of Galaxy objects). */
class Survey {
public:
    Survey() {}
    virtual~ Survey() {}

    /* Return the coordinate system for this survey. */
    virtual int GetCoordinateSystem() const = 0;

    /* Return a SeparationFunc object for this survey. */
    virtual SeparationFunc GetSeparationFunction() = 0;

    /* Return a callable SelectionFunc object representing the selection
     * function (expected number density of galaxies) for this survey. */
    virtual SelectionFunc GetSelectionFunction() = 0;

    /* Populate the vector 'gals' with the galaxies in the survey.  If the
     * supplied vector is non-empty the galaxies are added at the end, without
     * truncating it. */
    virtual void GetGalaxies(std::vector<Galaxy>& gals) = 0;
};

#include "cfg.h"

/* Load a concrete Survey instance based on the configuration options
 * specified in 'cfg'.  The configuration option "survey" must be set, and
 * must correspond to a fully implemented subclass of Survey. */
Survey* InitializeSurvey(Config cfg);

#endif // SURVEY_H
