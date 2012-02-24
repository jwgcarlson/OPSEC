#ifndef SEPARATIONFUNC_H
#define SEPARATIONFUNC_H

#include "opsec.h"

struct Point;
struct SeparationFuncImpl;

/* SeparationFunc
 *
 * Abstraction for quantifying the separation between two points in various
 * parameterizations. */
class SeparationFunc {
public:
    SeparationFunc(SeparationFuncImpl* impl = NULL);
    SeparationFunc(const SeparationFunc& sep);
    SeparationFunc& operator=(const SeparationFunc& sep);
    ~SeparationFunc();

    /* Return the comoving distance r between two points. */
    double r(const Point& p1, const Point& p2) const;

    /* Return the comoving distance r and the cosine mu of the angle between
     * the pair separation and the line-of-sight direction. */
    void rmu(const Point& p1, const Point& p2, double& r, double& mu) const;

    /* Return the comiving distance r between two points, as well as the
     * distances r1 and r2 of each point from the origin. */
    void rab(const Point& p1, const Point& p2, double& r, double& r1, double& r2) const;

private:
    SeparationFuncImpl* impl;
};

struct SeparationFuncImpl {
    /* Concrete implementations must implement the following methods */
    virtual double r(const Point& p1, const Point& p2) = 0;
    virtual void rmu(const Point& p1, const Point& p2, double& r, double& mu) = 0;
    virtual void rab(const Point& p1, const Point& p2, double& r, double& a, double& b) = 0;

protected:
    friend class SeparationFunc;
    int refcount;
    SeparationFuncImpl() { refcount = 0; }
    virtual ~SeparationFuncImpl() {}
};

#endif // SEPARATIONFUNC_H
