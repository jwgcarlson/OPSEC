#include <cstdio>
#include <cstring>

#include "Survey.h"


#include "BoxSurvey.h"
#include "SphericalSurvey.h"

Survey* InitializeSurvey(Config cfg) {
    const char* survey = cfg_get(cfg, "survey");
    Survey* s = NULL;
    Config subcfg = NULL;

#define HANDLE_SURVEY(NAME) \
    if(strcmp(survey, #NAME) == 0) { \
        subcfg = cfg_new_sub(cfg, #NAME ".", 1); \
        s = new NAME(subcfg); \
    }

    HANDLE_SURVEY(BoxSurvey)
    else HANDLE_SURVEY(SphericalSurvey)
    else {
        fprintf(stderr, "InitializeSurvey: unrecognized survey '%s'\n", survey);
        fflush(stderr);
    }
#undef HANDLE_SURVEY

    cfg_destroy(subcfg);
    return s;
}
