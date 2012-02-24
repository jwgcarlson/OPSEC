#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include "XiFunc.h"
#include "Cell.h"
#include "SeparationFunc.h"

XiFunc::XiFunc(XiFuncImpl* impl_) {
    impl = impl_;
    if(impl)
        ++impl->refcount;
}

XiFunc::XiFunc(const XiFunc& xi) {
    impl = xi.impl;
    if(impl)
        ++impl->refcount;
}

XiFunc& XiFunc::operator=(const XiFunc& xi) {
    if(impl != xi.impl) {
        if(impl && --impl->refcount <= 0)
            delete impl;
        impl = xi.impl;
        if(impl)
            ++impl->refcount;
    }
    return *this;
}

XiFunc::~XiFunc() {
    if(impl && --impl->refcount <= 0)
        delete impl;
}

double XiFunc::operator()(const Point& p1, const Point& p2, const SeparationFunc& sep) const {
#ifdef OPSEC_DEBUG
    if(impl == NULL) fprintf(stderr, "[XiFunc::operator()]: impl = NULL\n");
#endif
    return impl->xi(p1, p2, sep);
}
