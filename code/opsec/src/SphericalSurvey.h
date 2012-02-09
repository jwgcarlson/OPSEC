#ifndef SPHERICALSURVEY_H
#define SPHERICALSURVEY_H

#include <string>

#include "Survey.h"

/* Required configuration options:
 *   galfile
 *   maskfile
 *   radial
 */

class SphericalSurvey : public Survey {
public:
    SphericalSurvey(Config cfg);
    ~SphericalSurvey();

    SelectionFunc GetSelectionFunction();
    Galaxy* GetGalaxies(int* ngals);

    /* Compute the sky area coverage of the survey, weighted by the
     * completeness mask. */
    static double ComputeWeightedArea(long nside, float* mask);

    static void ComputeRadialNbar(int ngals, const Galaxy* gals, double A, int nr, const double* r, double* nbar);

protected:
    std::string galfile;
    std::string maskfile;
    std::string radial;
};

#endif // SPHERICALSURVEY_H
