#include <cstdio>
#include <cstring>

#include "Survey.h"

SelectionFunc::SelectionFunc(SelectionFuncImpl* impl_) {
    impl = impl_;
    if(impl)
        impl->refcount++;
}

SelectionFunc::SelectionFunc(const SelectionFunc& nbar) {
    impl = nbar.impl;
    if(impl)
        impl->refcount++;
}

SelectionFunc& SelectionFunc::operator=(const SelectionFunc& nbar) {
    if(impl != nbar.impl) {
        if(impl && --impl->refcount <= 0)
            delete impl;
        impl = nbar.impl;
        if(impl)
            impl->refcount++;
    }
    return *this;
}

SelectionFunc::~SelectionFunc() {
    if(impl && --impl->refcount <= 0)
        delete impl;
}

double SelectionFunc::operator()(double x, double y, double z) {
#ifdef DEBUG
    return impl ? impl->evaluate(x, y, z) : 0.0;
#else
    return impl->evaluate(x, y, z);
#endif
}


Survey::Survey() {
}

Survey::~Survey() {
}


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
