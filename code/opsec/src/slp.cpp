#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <algorithm>

#include "cscalapack.h"
#include "slp.h"

namespace slp {

Context::Context(int numrows, int numcols, char order, int sysctxt) {
    if(sysctxt == -1)
        Cblacs_get(0, 0, &ictxt);
    else
        ictxt = sysctxt;
    Cblacs_gridinit(&ictxt, (order == 'C') ? "Col" : "Row", numrows, numcols);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    if(myrow >= numrows || mycol >= numcols)
        myrow = mycol = -1;
}

Context::~Context() {
}

void Context::exit() {
    Cblacs_gridexit(ictxt);
    ictxt = nprow = npcol = myrow = mycol = -1;
}

int Context::numrows(int m, int mb, int rsrc) const {
    return numroc_(&m, &mb, &myrow, &rsrc, &nprow);
}

int Context::numcols(int n, int nb, int csrc) const {
    return numroc_(&n, &nb, &mycol, &csrc, &npcol);
}

Descriptor* Context::new_descriptor(int m, int n, int mb, int nb, int rsrc, int csrc, int lld) const {
    Descriptor* desc = new Descriptor();
    desc->dtype = 1;
    desc->ictxt = ictxt;
    desc->m = m;
    desc->n = n;
    desc->mb = mb;
    desc->nb = nb;
    desc->rsrc = rsrc;
    desc->csrc = csrc;
    desc->context = this;
    desc->mloc = numrows(m, mb, rsrc);
    desc->nloc = numcols(n, nb, csrc);
    desc->lld = std::max(lld, desc->mloc);
    return desc;
}

Descriptor::Descriptor() {
//    dtype = ictxt = m = n = mb = nb = rsrc = csrc = lld = mloc = nloc = -1;
//    context = NULL;
}

Descriptor::~Descriptor() {
}

int Descriptor::num_local_rows() const {
    return mloc;
}

int Descriptor::num_local_cols() const {
    return nloc;
}

int Descriptor::local_size() const {
    return lld*nloc;
}

int Descriptor::row_l2g(int iloc) const {
    int myrow = context->myrow;
    int nprow = context->nprow;
    return nprow*mb*(iloc/mb) + (iloc % mb) + ((nprow + myrow - rsrc) % nprow)*mb;
//    iloc += 1;
//    return indxl2g_(&iloc, &mb, &context->myrow, &rsrc, &context->nprow) - 1;
}

int Descriptor::col_l2g(int jloc) const {
    int mycol = context->mycol;
    int npcol = context->npcol;
    return npcol*nb*(jloc/nb) + (jloc % nb) + ((npcol + mycol - csrc) % npcol)*nb;
//    jloc += 1;
//    return indxl2g_(&jloc, &nb, &context->mycol, &csrc, &context->npcol) - 1;
}

void multiply(const Matrix<float>& A, const Matrix<float>& B, Matrix<float>& C,
        char transa, char transb, int m, int n, int k, float alpha, float beta,
        int ia, int ja, int ib, int jb, int ic, int jc)
{
    if(m <= 0)
        m = (transa == 'N') ? A.desc->m : A.desc->n;
    if(n <= 0)
        n = (transb == 'N') ? B.desc->n : B.desc->m;
    if(k <= 0)
        k = (transa == 'N') ? A.desc->n : A.desc->m;
    assert(C.desc->m >= m && C.desc->n >= n);
    ia += 1; ja += 1; ib += 1; jb += 1; ic += 1; jc += 1;       // Fortran indexing convention
    psgemm_(&transa, &transb, &m, &n, &k, &alpha, A.values, &ia, &ja, (int*) A.desc, B.values, &ib, &jb, (int*) B.desc, &beta, C.values, &ic, &jc, (int*) C.desc);
}

void multiply(const Matrix<double>& A, const Matrix<double>& B, Matrix<double>& C,
        char transa, char transb, int m, int n, int k, double alpha, double beta,
        int ia, int ja, int ib, int jb, int ic, int jc)
{
    if(m <= 0)
        m = (transa == 'N') ? A.desc->m : A.desc->n;
    if(n <= 0)
        n = (transb == 'N') ? B.desc->n : B.desc->m;
    if(k <= 0)
        k = (transa == 'N') ? A.desc->n : A.desc->m;
    assert(C.desc->m >= m && C.desc->n >= n);
    ia += 1; ja += 1; ib += 1; jb += 1; ic += 1; jc += 1;       // Fortran indexing convention
    pdgemm_(&transa, &transb, &m, &n, &k, &alpha, A.values, &ia, &ja, (int*) A.desc, B.values, &ib, &jb, (int*) B.desc, &beta, C.values, &ic, &jc, (int*) C.desc);
}

void redistribute(int m, int n, const Matrix<double>& A, int ia, int ja, Matrix<double>& B, int ib, int jb) {
    int ictxt;
    Cblacs_get(0, 0, &ictxt);
    Cpdgemr2d(m, n,
              (double*) A.values, ia + 1, ja + 1, (int*) A.desc, 
              B.values, ib + 1, jb + 1, (int*) B.desc, 
              ictxt);
}

void redistribute(int m, int n, const Matrix<double>& A, int ia, int ja, Matrix<double>& B, int ib, int jb, const Context& gcontext) {
    Cpdgemr2d(m, n,
              (double*) A.values, ia + 1, ja + 1, (int*) A.desc,
              B.values, ib + 1, jb + 1, (int*) B.desc, 
              gcontext.ictxt);
}

void gsum2d(const Context* context, const char* scope, const char* top,
        int m, int n, float* a, int lda, int rdest, int cdest)
{
    Csgsum2d(context->ictxt, (char*) scope, (char*) top, m, n, (char*) a, lda, rdest, cdest);
}

void gsum2d(const Context* context, const char* scope, const char* top,
        int m, int n, double* a, int lda, int rdest, int cdest)
{
    Cdgsum2d(context->ictxt, (char*) scope, (char*) top, m, n, (char*) a, lda, rdest, cdest);
}

} // namespace slp




#if 0
int myindxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs) {
    iloc += 1;
    return indxl2g_(&iloc, &nb, &iproc, &isrcproc, &nprocs) - 1;
}

void mypdgemv(const char* transa, int m, int n,
              double alpha,
              const double* a, int ia, int ja, const ArrayDesc& desca,
              const double* x, int ix, int jx, const ArrayDesc& descx, int incx,
              double beta,
              double* y, int iy, int jy, const ArrayDesc& descy, int incy)
{
    ia++; ja++; ix++; jx++; iy++; jy++;         // convert to Fortran indexing
    pdgemv_((char*) transa, &m, &n,
            &alpha,
            (double*) a, &ia, &ja, (int*) &desca,
            (double*) x, &ix, &jx, (int*) &descx, &incx,
            &beta,
            y, &iy, &jy, (int*) &descy, &incy);
}

void mypdgemm(const char* transa, const char* transb, int m, int n, int k,
              double alpha,
              const double* a, int ia, int ja, const ArrayDesc& desca,
              const double* b, int ib, int jb, const ArrayDesc& descb,
              double beta,
              double* c, int ic, int jc, const ArrayDesc& descc)
{
    ia++; ja++; ib++; jb++; ic++; jc++;         // convert to Fortran indexing
    pdgemm_((char*) transa, (char*) transb, &m, &n, &k,
            &alpha,
            (double*) a, &ia, &ja, (int*) &desca,
            (double*) b, &ib, &jb, (int*) &descb,
            &beta,
            (double*) c, &ic, &jc, (int*) &descc);
}

void mypdgemr2d(int m, int n,
                const double* a, int ia, int ja, const ArrayDesc& desca,
                double* b, int ib, int jb, const ArrayDesc& descb,
                const Context& gcontext)
{
    Cpdgemr2d(m, n,
              (double*) a, ia + 1, ja + 1, (int*) &desca, 
              b, ib + 1, jb + 1, (int*) &descb, 
              gcontext.ictxt);
}
#endif
