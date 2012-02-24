#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include "SelectionFunc.h"
#include "Cell.h"

SelectionFunc::SelectionFunc(SelectionFuncImpl* impl_) {
    impl = impl_;
    if(impl)
        ++impl->refcount;
}

SelectionFunc::SelectionFunc(const SelectionFunc& nbar) {
    impl = nbar.impl;
    if(impl)
        ++impl->refcount;
}

SelectionFunc& SelectionFunc::operator=(const SelectionFunc& nbar) {
    if(impl != nbar.impl) {
        if(impl && --impl->refcount <= 0)
            delete impl;
        impl = nbar.impl;
        if(impl)
            ++impl->refcount;
    }
    return *this;
}

SelectionFunc::~SelectionFunc() {
    if(impl && --impl->refcount <= 0)
        delete impl;
}

double SelectionFunc::operator()(double x, double y, double z) {
#ifdef DEBUG
   if(impl == NULL) fprintf(stderr, "[SelectionFunc::operator()]: impl = NULL\n"); 
#endif
    return impl->evaluate(x, y, z);
}
