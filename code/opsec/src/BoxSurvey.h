#ifndef BOXSURVEY_H
#define BOXSURVEY_H

#include <string>

#include "Survey.h"

/* Required configuration options:
 *   galfile
 *   nbar
 *   L
 */

class BoxSurvey : public Survey {
public:
    BoxSurvey(const char* galfile, double nbar, double L);
    BoxSurvey(Config cfg);
    ~BoxSurvey();

    /* Survey interface */
    int GetCoordinateSystem() const;
    SeparationFunc GetSeparationFunction();
    SelectionFunc GetSelectionFunction();
    void GetGalaxies(std::vector<Galaxy>& gals);

protected:
    std::string galfile;
    double nbar;        // mean galaxy number density
    double L;           // box size
};

#endif // BOXSURVEY_H
