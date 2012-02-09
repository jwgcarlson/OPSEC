#ifndef BOXSURVEY_H
#define BOXSURVEY_H

#include <string>

#include "Survey.h"

/* Required configuration options:
 *   galfile
 *   nbar
 */

class BoxSurvey : public Survey {
public:
    BoxSurvey(Config cfg);
    ~BoxSurvey();

    SelectionFunc GetSelectionFunction();
    Galaxy* GetGalaxies(int* ngals);

protected:
    std::string galfile;
    double nbar;
};

#endif // BOXSURVEY_H
