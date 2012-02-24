#ifndef SELECTIONFUNC_H
#define SELECTIONFUNC_H

#include "opsec.h"

class Point;
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
    double operator()(double x1, double x2, double x3);
    double operator()(const Point& p);

private:
    SelectionFuncImpl* impl;
};

struct SelectionFuncImpl {
    /* Concrete implementations must define this method */
    virtual double evaluate(double x1, double x2, double x3) = 0;

protected:
    friend class SelectionFunc;
    int refcount;
    SelectionFuncImpl() { refcount = 0; }
    virtual ~SelectionFuncImpl() {}
};

#endif // SELECTIONFUNC_H
