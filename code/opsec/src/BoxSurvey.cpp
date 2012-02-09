#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "BoxSurvey.h"
#include "Spline.h"
#include "abn.h"


/* Constant selection function. */
struct ConstantSelectionFuncImpl : public SelectionFuncImpl {
    double const_nbar;          // constant value for $\bar{n}(\vec{r})$

    ConstantSelectionFuncImpl(double const_nbar_) {
        const_nbar = const_nbar_;
        assert(const_nbar > 0);
    }

    ~ConstantSelectionFuncImpl() {}

    double evaluate(double, double, double) {
        return const_nbar;
    }
};


BoxSurvey::BoxSurvey(Config cfg) {
    galfile = cfg_get(cfg, "galfile");
    nbar = cfg_get_double(cfg, "nbar");
}

BoxSurvey::~BoxSurvey() {
}

SelectionFunc BoxSurvey::GetSelectionFunction() {
    return new ConstantSelectionFuncImpl(nbar);
}

Galaxy* BoxSurvey::GetGalaxies(int* ngals) {
    FILE* fgals = fopen(galfile.c_str(), "r");
    if(fgals == NULL) {
        fprintf(stderr, "BoxSurvey: could not open file '%s'\n", galfile.c_str());
        return NULL;
    }

    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    if(abn_read_header(fgals, &n, &size, &endian, fmt, NULL) != 0
       || size != sizeof(Galaxy)
       || strcmp(fmt, "4f") != 0)
    {
        fprintf(stderr, "BoxSurvey: error reading galaxies from '%s'\n", galfile.c_str());
        fclose(fgals);
        return NULL;
    }

    Galaxy* gals = (Galaxy*) malloc(n*sizeof(Galaxy));
    if(!gals) {
        fprintf(stderr, "BoxSurvey: could not allocate memory for %zd galaxies\n", n);
        fclose(fgals);
        return NULL;
    }

    size_t nread = fread((void*) gals, size, n, fgals);
    if(nread != n) {
        fprintf(stderr, "BoxSurvey: expecting %zd galaxies from '%s', got %zd\n", n, galfile.c_str(), nread);
        fclose(fgals);
        return NULL;
    }

    fclose(fgals);
    if(ngals)
        *ngals = (int)n;
    return gals;
}
