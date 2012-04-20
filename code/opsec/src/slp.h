#ifndef SLP_H
#define SLP_H

#include <cassert>

#include "cscalapack.h"

namespace slp {

struct Context;
struct Descriptor;

/* BLACS context */
struct Context {
    int ictxt;      // BLACS context ID
    int nprow;      // number of rows in process grid
    int npcol;      // number of columns in process grid
    int myrow;      // row number of this process (0 <= myrow < nprow)
    int mycol;      // column number of this process (0 <= mycol < npcol)

public:
    /* Establish a nrows-by-ncols process grid, using either row-major (order =
     * 'R') or column-major (order = 'C') ordering.  If not equal to -1,
     * sysctxt indicates the system context to be used in creating the new
     * BLACS context; if not specified, the default system context is used. */
    Context(int nrows, int ncols, char order = 'R', int sysctxt = -1);

    /* Does nothing.  Context must be deallocated manually with exit(). */
    ~Context();

    /* Exit context and clean up associated resources. */
    void exit();

    /* Return the number of local rows or columns that would belong to the
     * current process, for a m-by-n matrix distributed as mb-by-nb blocks. */
    int numrows(int m, int mb, int rsrc = 0) const;
    int numcols(int n, int nb, int csrc = 0) const;

    /* Create a new array descriptor for m-by-n matrices broken down into
     * mb-by-nb blocks.  (rsrc,csrc) denote the coordinates of the process that
     * should hold the upper-left block of the matrix.  lld is the leading
     * dimension of the local array that will hold the matrix elements; 0
     * indicates that a default value should be chosen. */
    Descriptor* new_descriptor(int m, int n, int mb, int nb,
                               int rsrc = 0, int csrc = 0, int lld = 0) const;
};

/* BLACS array descriptor */
struct Descriptor {
    int dtype;  // = 1 for dense matrices
    int ictxt;  // BLACS context handle
    int m;      // number of rows
    int n;      // number of columns
    int mb;     // row blocking factor
    int nb;     // column blocking factor
    int rsrc;   // = 0
    int csrc;   // = 0
    int lld;    // leading dimension of local array

    int mloc;   // number of local rows
    int nloc;   // number of local columns

public:
    ~Descriptor();

    /* Return int array for use in Fortran routines expecting a descriptor */
    int* intarray() const;

    int num_local_rows() const;
    int num_local_cols() const;

    /* Number of elements to be allocated locally, locsize = lld*nloc */
    int local_size() const;

    int row_l2g(int iloc) const;
    int col_l2g(int jloc) const;

private:
    /* Descriptors should only be allocated by a Context. */
    friend class Context;
    Descriptor();

    const Context* context;     // BLACS context object (TODO: make this a shared pointer)
};

struct LocalMatrixIndex {
    int i, j;
    int iloc, jloc;

    LocalMatrixIndex(const Descriptor* desc_, int kloc = 0) {
        desc = desc_;
        set(kloc);
    }

    /* To allow initialization from a Matrix object */
    template<class T>
    LocalMatrixIndex(const T& x, int kloc = 0) {
        desc = x.desc;
        set(kloc);
    }

    operator int() const { return iloc + jloc*desc->lld; }

    bool done() const {
        return (kloc == desc->mloc*desc->nloc);
    }

    void set(int kloc_) {
        assert(0 <= kloc_ && kloc_ <= desc->mloc*desc->nloc); // values outside this range are most likely an error
        kloc = kloc_;
        iloc = kloc % desc->mloc;
        jloc = kloc/desc->mloc;
        i = desc->row_l2g(iloc);
        j = desc->col_l2g(jloc);
    }

    LocalMatrixIndex& operator=(int kloc_) { set(kloc_);     return *this; }
    LocalMatrixIndex& operator+=(int dk)   { set(kloc + dk); return *this; }
    LocalMatrixIndex& operator-=(int dk)   { set(kloc - dk); return *this; }
    LocalMatrixIndex& operator++()         { set(kloc + 1);  return *this; }
    LocalMatrixIndex& operator--()         { set(kloc - 1);  return *this; }

private:
    int kloc;           // 0 <= kloc < mloc*nloc, iloc = kloc % mloc, jloc = kloc/mloc
    const Descriptor* desc;
};

template<class T>
struct Matrix {
    const Descriptor* desc;
    T* values;

public:
    Matrix() {
        desc = NULL;
        values = NULL;
    }

    Matrix(const Descriptor* desc_, T* values_) : desc(desc_), values(values_) {
        assert(desc != NULL && values != NULL);
    }

    ~Matrix() {}

    /* Flat local accessors (for vectors or flattened matrices) */
    T& operator[](int k) { return values[k]; }
    const T& operator[](int k) const { return values[k]; }

    /* 2-D local accessors */
    T& operator()(int iloc, int jloc) { return values[iloc + desc->lld*jloc]; }
    const T& operator()(int iloc, int jloc) const { return values[iloc + desc->lld*jloc]; }

#if 0
    /* Iterator class */
    template<class PT>
    struct MatrixIterator {
        int &i, &j;
        int &iloc, &jloc;

        MatrixIterator(const Descriptor* desc, PT values_, int kset = 0)
            : k(desc, kset), values(values_), i(k.i), j(k.j), iloc(k.iloc), jloc(k.jloc)
        {}

        PT operator*() { return &values[k]; }

        MatrixIterator& operator++()       {     ++k; return *this; }
        MatrixIterator& operator--()       {     --k; return *this; }
        MatrixIterator& operator+=(int dk) { k += dk; return *this; }
        MatrixIterator& operator-=(int dk) { k -= dk; return *this; }
        bool operator<(PT x) { return (this->operator() < x); }
        bool operator<=(PT x) { return (this->operator() <= x); }
        bool operator==(PT x) { return (this->operator() == x); }
        bool operator!=(PT x) { return (this->operator() != x); }
        bool operator>=(PT x) { return (this->operator() >= x); }
        bool operator>(PT x) { return (this->operator() > x); }

    private:
        LocalMatrixIndex k;
        PT values;
    };

    typedef MatrixIterator<T*> iterator;
    typedef MatrixIterator<const T*> const_iterator;

    iterator begin() { return iterator(desc, values, 0); }
    const_iterator begin() const { return const_iterator(desc, values, 0); }
    iterator end() { return iterator(desc, values, mloc*nloc); }
    const_iterator end() const { return const_iterator(desc, values, mloc*nloc); }
#endif
};

/* Matrix-matrix multiplication.  These are thin wrappers around psgemm() and
 * pdgemm(). */
void multiply(const Matrix<float>& A, const Matrix<float>& B, Matrix<float>& C,
              char transa = 'N', char transb = 'N',
              int m = 0, int n = 0, int k = 0,
              float alpha = 1.0f, float beta = 0.0f,
              int ia = 0, int ja = 0, int ib = 0, int jb = 0, int ic = 0, int jc = 0);
void multiply(const Matrix<double>& A, const Matrix<double>& B, Matrix<double>& C,
              char transa = 'N', char transb = 'N',
              int m = 0, int n = 0, int k = 0,
              double alpha = 1.0, double beta = 0.0,
              int ia = 0, int ja = 0, int ib = 0, int jb = 0, int ic = 0, int jc = 0);

/* Matrix redistribution */
void redistribute(int m, int n,
                  const Matrix<double>& A, int ia, int ja,
                  Matrix<double>& B, int ib, int jb);
void redistribute(int m, int n,
                  const Matrix<double>& A, int ia, int ja,
                  Matrix<double>& B, int ib, int jb,
                  const Context& gcontext);

/* Element-wise summation */
/* TODO: wrap this in a more intuitive interface? */
void gsum2d(const Context* context, const char* scope, const char* top,
            int m, int n, float* a, int lda,
            int rdest = -1, int cdest = -1);
void gsum2d(const Context* context, const char* scope, const char* top,
            int m, int n, double* a, int lda,
            int rdest = -1, int cdest = -1);





#if 0
Context create_context(int numrows, int numcols);

int mynumroc(int n, int nb, int iproc, int isrcproc, int nprocs);

ArrayDesc create_descriptor(const Context& context,
                            int m, int n, int mb, int nb,
                            int rsrc = 0, int csrc = 0, int lld = 0);

int myindxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs);

void mypdgemv(const char* transa, int m, int n,
              double alpha,
              const double* a, int ia, int ja, const ArrayDesc& desca,
              const double* x, int ix, int jx, const ArrayDesc& descx, int incx,
              double beta,
              double* y, int iy, int jy, const ArrayDesc& descy, int incy);

void mypdgemm(const char* transa, const char* transb, int m, int n, int k,
              double alpha,
              const double* a, int ia, int ja, const ArrayDesc& desca,
              const double* b, int ib, int jb, const ArrayDesc& descb,
              double beta,
              double* c, int ic, int jc, const ArrayDesc& descc);

void mypdgemr2d(int m, int n,
                const double* a, int ia, int ja, const ArrayDesc& desca,
                double* b, int ib, int jb, const ArrayDesc& descb,
                const Context& gcontext);
#endif

} // namespace slp

#endif // SLP_H
