#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include "SeparationFunc.h"
#include "Cell.h"

SeparationFunc::SeparationFunc(SeparationFuncImpl* impl_) {
    impl = impl_;
    if(impl)
        ++impl->refcount;
}

SeparationFunc::SeparationFunc(const SeparationFunc& sep) {
    impl = sep.impl;
    if(impl)
        ++impl->refcount;
}

SeparationFunc& SeparationFunc::operator=(const SeparationFunc& sep) {
    if(impl != sep.impl) {
        if(impl && --impl->refcount <= 0)
            delete impl;
        impl = sep.impl;
        if(impl)
            ++impl->refcount;
    }
    return *this;
}

SeparationFunc::~SeparationFunc() {
    if(impl && --impl->refcount <= 0)
        delete impl;
}

double SeparationFunc::r(const Point& p1, const Point& p2) const {
#ifdef OPSEC_DEBUG
    if(impl == NULL) fprintf(stderr, "[SeparationFunc::r]: impl = NULL\n");
#endif
    return impl->r(p1, p2);
}

void SeparationFunc::rmu(const Point& p1, const Point& p2, double& r, double& mu) const {
#ifdef OPSEC_DEBUG
    if(impl == NULL) fprintf(stderr, "[SeparationFunc::rmu]: impl = NULL\n");
#endif
    impl->rmu(p1, p2, r, mu);
}

void SeparationFunc::rab(const Point& p1, const Point& p2, double& r, double& a, double& b) const {
#ifdef OPSEC_DEBUG
    if(impl == NULL) fprintf(stderr, "[SeparationFunc::rab]: impl = NULL\n");
#endif
    impl->rab(p1, p2, r, a, b);
}
