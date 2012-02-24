#ifndef XIFUNC_H
#define XIFUNC_H

#include <cassert>

#include "opsec.h"

class Point;
class SeparationFunc;
struct XiFuncImpl;

/* A XiFunc object represents a 2-point function, i.e. a real-valued function
 * that depends on the location of two points p1 and p2.  Its primary use is in
 * computing the components of the signal matrix.
 *
 * The XiFunc class is a thin wrapper around a XiFuncImpl object, which
 * actually performs the task of computing the 2-point signal.  This pattern
 * allows for abstraction (XiFuncImpl is an abstract class) while avoiding
 * having to explicitly pass around pointers or references (and the resulting
 * confusion with regard to object permanence and ownership).  Overall I just
 * think it makes the code look a little cleaner. */
class XiFunc {
public:
    XiFunc(XiFuncImpl* impl = NULL);
    XiFunc(const XiFunc& xi);
    XiFunc& operator=(const XiFunc& xi);
    ~XiFunc();

    double operator()(const Point& p1, const Point& p2, const SeparationFunc& sep) const;

private:
    XiFuncImpl* impl;
};

struct XiFuncImpl {
    /* Concrete implementations must implement the following method */
    virtual double xi(const Point& p1, const Point& p2, const SeparationFunc& sep) = 0;

protected:
    friend class XiFunc;
    int refcount;
    XiFuncImpl() { refcount = 0; }
    virtual ~XiFuncImpl() { assert(refcount == 0); }
};

#endif // XIFUNC_H
