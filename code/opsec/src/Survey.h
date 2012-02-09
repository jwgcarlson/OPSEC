#ifndef SURVEY_H
#define SURVEY_H

#include "cfg.h"

/* Data structure representing a galaxy in either Cartesian or spherical coordinates. */
struct Galaxy {
    union { float x; float r; };
    union { float y; float mu; };
    union { float z; float phi; };
    float w;    // galaxy weight
};

struct SelectionFuncImpl;

/* SelectionFunc
 * 
 * A callable object that evaluates the survey's selection function.  All OPSEC
 * routines that require a selection function take an object of this type.  The
 * details of how the selection function is defined are encapsulated in the
 * SelectionFuncImpl class. */
class SelectionFunc {
public:
    SelectionFunc(SelectionFuncImpl* impl = NULL);
    SelectionFunc(const SelectionFunc& nbar);
    SelectionFunc& operator=(const SelectionFunc& nbar);
    ~SelectionFunc();

    /* Evaluate the selection function at the specified point [either (x,y,z) 
     * for Cartesian or (r,mu,phi) for spherical coordinates] */
    double operator()(double, double, double);

private:
    SelectionFuncImpl* impl;
};

/* Implementation of a selection function. */
struct SelectionFuncImpl {
    virtual double evaluate(double, double, double) = 0;

    int refcount;
    SelectionFuncImpl() { refcount = 0; }
    virtual ~SelectionFuncImpl() {}
};


/* Survey
 * 
 * Abstract represention of a galaxy survey, including the selection function
 * (expected number density of galaxies) and the observed galaxy positions. */
class Survey {
public:
    Survey();
    virtual~ Survey();

    /* Return a callable SelectionFunc object representing the selection
     * function (expected number density of galaxies) for this survey. */
    virtual SelectionFunc GetSelectionFunction() = 0;

    /* Return a pointer to an array of galaxies in the survey.  The array is
     * allocated with malloc(), so the caller must free() the pointer when done
     * with the galaxies.  'ngals' is set to the size of the array. */
    virtual Galaxy* GetGalaxies(int* ngals) = 0;
};

/* Load a concrete Survey instance based on the configuration options
 * specified in 'cfg'.  The configuration option "survey" must be set, and
 * must correspond to a fully implemented subclass of Survey. */
Survey* InitializeSurvey(Config cfg);

#endif // SURVEY_H
