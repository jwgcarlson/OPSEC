#ifndef MYMATRIX_HPP
#define MYMATRIX_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <Teuchos_BLAS.hpp>
#include <Teuchos_LAPACK.hpp>

/* Simple classes for treating regular arrays of numbers as vectors or matrices. */

template<class T>
static inline T mysum(int n, T* x, int incx) {
    T sum = 0;
    for(int i = 0; i < n; i++) {
        sum += (*x);
        x += incx;
    }
    return sum;
}

/* Convenience class for keeping track of indices for packed matrices.
 * Currently it only handles upper triangular (column-major) matrices, but it
 * can be trivially modified to support other cases. */
struct PackedMatrixIndex {
    explicit PackedMatrixIndex(int n_) {
        n = n_;
        k = i = j = 0;
    }

    PackedMatrixIndex(int n_, int k_) {
        n = n_;
        set(k_);
    }

    void set(int knew) {
        k = knew;
        i = j = 0;
        while(j*(j+1)/2 < k)
            ++j;
        while(j*(j+1)/2 + i < k)
            ++i;
    }

    void update(int knew) {
        while(k < knew) {
            /* Count up to knew */
            ++i;
            if(i > j) {
                ++j;
                i = 0;
            }
            ++k;
        }
        while(k > knew) {
            /* Count down to knew */
            --i;
            if(i < 0) {
                --j;
                i = j;
            }
            --k;
        }
    }

    PackedMatrixIndex& operator=(int knew) { set(knew); return *this; }
    PackedMatrixIndex& operator+=(int dk) { update(k + dk); return *this; }
    PackedMatrixIndex& operator-=(int dk) { update(k - dk); return *this; }
    PackedMatrixIndex& operator++() { update(k + 1); return *this; }
    PackedMatrixIndex& operator--() { update(k - 1); return *this; }
    operator int() { return k; }
//    bool operator<(int x) { return (k < x); }
//    bool operator<=(int x) { return (k <= x); }
//    bool operator==(int x) { return (k == x); }
//    bool operator>=(int x) { return (k >= x); }
//    bool operator>(int x) { return (k > x); }

    int n;      // size of matrix (n-by-n)
    int k;      // index into packed array of matrix elements [0 <= k < n*(n+1)/2]
    int i, j;   // matrix indices corresponding to packed index k
};

/* Matrix storage order */
enum {
    ColumnMajor = 0,
    RowMajor = 1
};

/* MyVector
 *
 * Wrap an array of ScalarType values as a vector. */
template<class ScalarType>
struct MyVector : public Teuchos::BLAS<int,ScalarType> {
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;


    /***** Member variables **************************************************/

    int n;              // length of vector
    int inc;            // increment between adjacent elements
    ScalarType* values;


    /****** Constructors *****************************************************/

    MyVector() {
        n = inc = 0;
        values = NULL;
    }

    MyVector(int n_, ScalarType* values_, int inc_ = 1) {
        n = n_;
        inc = inc_;
        values = values_;
    }

    MyVector(const MyVector& y) {
        n = y.n;
        inc = y.inc;
        values = y.values;
    }

    ~MyVector() {
    }

    MyVector& operator=(const MyVector& y) {
        n = y.n;
        inc = y.inc;
        values = y.values;
        return *this;
    }

    /* Copy this vector into the vector y. */
    void copy(MyVector& y) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && y.values != NULL && n == y.n);
#endif
        COPY(n, values, inc, y.values, y.inc);
    }

    /* Copy this vector into a new vector.  The array newvalues must be able
     * to hold at least n elements, and must not overlap with this->values. */
    MyVector newcopy(ScalarType* newvalues) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && newvalues != NULL);
#endif
        COPY(n, values, inc, newvalues, 1);
        return MyVector(n, newvalues, 1);
    }


    MagnitudeType norm() const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && n > 0);
#endif
        return NRM2(n, values, inc);
    }

    void scale(ScalarType alpha) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && n > 0);
#endif
        SCAL(n, alpha, values, inc);
    }

    /***** Accessors **********************************************************/

    operator ScalarType*() {
        return values;
    }

    operator const ScalarType*() const {
        return values;
    }

    ScalarType& element(int i) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && 0 <= i && i < n);
#endif
        return values[i*inc];
    }

    ScalarType& operator[](int i) {
        return element(i);
    }

    const ScalarType& operator[](int i) const {
        return element(i);
    }

};

template<class ScalarType>
struct MyMatrix : public Teuchos::BLAS<int,ScalarType> {

    /***** Member variables **************************************************/

    int m;              // number of rows
    int n;              // number of columns
    int storage;        // storage order, either ColumnMajor or RowMajor
    int stride;         // column or row stride, depending on storage order
    ScalarType* values;


    /***** Constructors *******************************************************/

    /* Default constructor.  No methods may be called on a default-constructed
     * MyMatrix (they will probably cause a segfault). */
    MyMatrix() {
        m = n = stride = 0;
        storage = -1;
        values = NULL;
    }

    /* Set stride to 0 to indicate dense storage, i.e.
     *   stride = m if storage == ColumnMajor,
     *   stride = n if storage == RowMajor. */
    MyMatrix(int m_, int n_, ScalarType* values_, int storage_ = RowMajor, int stride_ = 0) {
        m = m_;
        n = n_;
        storage = storage_;
        stride = stride_;
        if(storage == ColumnMajor && stride < m)
            stride = m;
        else if(storage == RowMajor && stride < n)
            stride = n;
        values = values_;
    }

    /* Copy constructor. */
    MyMatrix(const MyMatrix& B) {
        m = B.m;
        n = B.n;
        storage = B.storage;
        stride = B.stride;
        values = B.values;
    }

    ~MyMatrix() {
    }

    MyMatrix& operator=(const MyMatrix& B) {
        m = B.m;
        n = B.n;
        storage = B.storage;
        stride = B.stride;
        values = B.values;
        return *this;
    }


    /* Amount of padding at the end of each row or column. */
    int padding() const {
        if(storage == ColumnMajor)
            return stride - m;
        else if(storage == RowMajor)
            return stride - n;
    }

    /* Create a view into a m1-by-n1 sub-matrix of this matrix. */
    MyMatrix block(int m1, int n1, int i1, int j1) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && m1 > 0 && n1 > 0 && 0 <= i1 && i1 + m1 <= m && 0 <= j1 && j1 + n1 <= n);
#endif
        return MyMatrix(m1, n1, &element(i1,j1), storage, stride);
    }

    /* Copy this matrix into the matrix B. */
    void copy(MyMatrix& B) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && B.values != NULL && A.m == B.m && A.n == B.n);
#endif
        if(padding() == 0 && B.padding() == 0) {
            /* Dense matrix storage, memcpy the entire array */
            memcpy(B.values, values, m*n*sizeof(ScalarType));
        }
        else if(storage == ColumnMajor) {
            /* Copy column by column */
            for(int j = 0; j < n; j++) {
                MyVector<ScalarType> a = this->column(j);
                MyVector<ScalarType> b = B.column(j);
                COPY(m, a, a.inc, b, b.inc);
            }
        }
        else if(storage == RowMajor) {
            /* Copy row by row */
            for(int i = 0; i < m; i++) {
                MyVector<ScalarType> a = this->row(i);
                MyVector<ScalarType> b = B.row(i);
                COPY(n, a, a.inc, b, b.inc);
            }
        }
    }

    /* Copy values into a new matrix.  The array newvalues must be able to hold
     * at least m*n elements, and must not overlap with this->values.  The new
     * matrix has dense storage, i.e. there is no padding after each row or
     * column. */
    MyMatrix newcopy(ScalarType* newvalues) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && newvalues != NULL);
#endif
        if((storage == ColumnMajor && stride == m) || (storage == RowMajor && stride == n)) {
            /* Dense matrix storage, memcpy the entire array */
            memcpy(newvalues, values, m*n*sizeof(ScalarType));
        }
        else if(storage == ColumnMajor) {
            /* Copy column by column */
            for(int j = 0; j < n; j++)
                memcpy(&newvalues[j*m], this->column(j), m*sizeof(ScalarType));
        }
        else {
            /* Copy row by row */
            for(int i = 0; i < m; i++)
                memcpy(&newvalues[i*n], this->row(i), n*sizeof(ScalarType));
        }
        int newstride = (storage == ColumnMajor) ? m : n;
        return MyMatrix(m, n, newvalues, storage, newstride);
    }

    /* Return matrix transpose. */
    MyMatrix transpose() {
        int newstorage = (storage == ColumnMajor) ? RowMajor : ColumnMajor;
        return MyMatrix(m, n, values, newstorage, stride);
    }

    /* Return trace of a square matrix. */
    ScalarType trace() {
#ifdef OPSEC_DEBUG
        assert(values != NULL && m == n);
#endif
        return mysum(n, values, stride+1);
    }

    /***** Accessors **********************************************************/

    operator ScalarType*() {
        return values;
    }

    operator const ScalarType*() const {
        return values;
    }

    /* Access element (i,j) of this matrix. */
    ScalarType& element(int i, int j) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && 0 <= i && i < m && 0 <= j && j < n);
#endif
        return (storage == ColumnMajor) ? values[i + j*stride]
                                      : values[i*stride + j];
    }

    const ScalarType& element(int i, int j) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && 0 <= i && i < m && 0 <= j && j < n);
#endif
        return (storage == ColumnMajor) ? values[i + j*stride]
                                      : values[i*stride + j];
    }

    ScalarType& operator()(int i, int j) {
        return element(i,j);
    }

    const ScalarType& operator()(int i, int j) const {
        return element(i,j);
    }

    MyVector<ScalarType> column(int j) const {
        if(storage == ColumnMajor)
            return MyVector<ScalarType>(m, &values[j*stride], 1);
        else
            return MyVector<ScalarType>(m, &values[j], stride);
    }

    MyVector<ScalarType> row(int i) const {
        if(storage == ColumnMajor)
            return MyVector<ScalarType>(n, &values[i], stride);
        else
            return MyVector<ScalarType>(n, &values[i*stride], 1);
    }


    /***** Input/output ******************************************************/

    void print(std::ostream& os = std::cout) const {
#ifdef OPSEC_DEBUG
        assert(values != NULL);
#endif
        for(int i = 0; i < m; i++) {
            os << ((i == 0) ? "[" : " ") << "[ ";
            for(int j = 0; j < n; j++)
                os << element(i,j) << " ";
            os << "]" << ((i == m-1) ? "]" : "") << std::endl;
        }
    }
};


/* Vector-vector addition */
template<class ScalarType>
void axpy(ScalarType alpha, const MyVector<ScalarType>& x, MyVector<ScalarType>& y) {
#ifdef OPSEC_DEBUG
    assert(x.values != NULL && y.values != NULL && x.n == y.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    blas.AXPY(x.n, alpha, x, x.inc, y, y.inc);
}

/* Vector inner product */
template<class ScalarType>
ScalarType dot(const MyVector<ScalarType>& x, const MyVector<ScalarType>& y) {
#ifdef OPSEC_DEBUG
    assert(x.values != NULL && y.values != NULL && x.n == y.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    return blas.DOT(x.n, x, x.inc, y, y.inc);
}

/* Matrix-vector product */
template<class ScalarType>
void multiply(const MyMatrix<ScalarType>& A, const MyVector<ScalarType>& x, MyVector<ScalarType>& y) {
#ifdef OPSEC_DEBUG
    assert(A.values != NULL && x.values != NULL && y.values != NULL && A.n == x.n && A.m == y.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::ETransp trans = (A.storage == ColumnMajor) ? Teuchos::NO_TRANS : Teuchos::TRANS;
    blas.GEMV(trans, A.m, A.n, 1, A, A.stride, x, x.inc, 0, y, y.inc);
}

/* Matrix-vector product when the matrix A is square diagonal, with $A_{ii} =
 * a_i$.  It is okay if y = x. */
template<class ScalarType>
void multiply_diagonal(const MyVector<ScalarType>& a, const MyVector<ScalarType>& x, MyVector<ScalarType>& y) {
#ifdef OPSEC_DEBUG
    assert(a.values != NULL && x.values != NULL && y.values != NULL && a.n == x.n && a.n == y.n);
#endif
    int n = a.n;
    #pragma omp parallel for
    for(int i = 0; i < n; i++)
        y[i] = a[i]*x[i];
}

/* Matrix-matrix product */
template<class ScalarType>
void multiply(const MyMatrix<ScalarType>& A, const MyMatrix<ScalarType>& B, MyMatrix<ScalarType>& C) {
#ifdef OPSEC_DEBUG
    assert(A.values != NULL && B.values != NULL && C.values != NULL && C.m == A.m && A.n == B.m && B.n == C.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::ETransp transa = (A.storage == ColumnMajor) ? Teuchos::NO_TRANS : Teuchos::TRANS;
    Teuchos::ETransp transb = (A.storage == ColumnMajor) ? Teuchos::NO_TRANS : Teuchos::TRANS;
    blas.GEMM(transa, transb, C.m, C.n, A.n, 1, A, A.stride, B, B.stride, 0, C, C.stride);
}

/* Matrix-matrix product when the matrix A is square diagonal, with $A_{ii} =
 * a_i$.  It is okay if B = C. */
template<class ScalarType>
void multiply_diagonal(const MyVector<ScalarType>& a, const MyMatrix<ScalarType>& B, MyMatrix<ScalarType>& C) {
#ifdef OPSEC_DEBUG
    assert(a.values != NULL && B.values != NULL && C.values != NULL && C.m == a.n && a.n == B.m && B.n == C.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    int m = B.m;        // number of rows
    int n = B.n;        // number of columns
    for(int i = 0; i < m; i++) {
        MyVector<ScalarType> b = B.row(i);
        MyVector<ScalarType> c = C.row(i);
        blas.COPY(n, b, b.inc, c, c.inc);
        c.scale(a[i]);
    }
}

/* Compute the principal square root of the symmetric, positive-definite matrix
 * A, storing the result in B.  It is okay if B points to the same memory as A. */
template<class ScalarType>
void sqrtm(const MyMatrix<ScalarType>& A, MyMatrix<ScalarType>& B) {
#ifdef OPSEC_DEBUG
    assert(A.values != NULL && A.m == A.n);
    assert(B.values != NULL && B.m == A.m && B.n == A.n);
#endif

    Teuchos::LAPACK<int,ScalarType> lapack;
    int n = A.n;
    int blocksize = lapack.ILAENV(1, "DSYTRD", "", n);
    int lwork = n*(blocksize+2);        // optimal work size
    char jobz = 'V', uplo = 'U';
    int info;

    /* Allocate temporary storage */
    ScalarType* w = (ScalarType*) malloc(n*sizeof(ScalarType)); // eigenvalues
    ScalarType* work = (ScalarType*) malloc(lwork*sizeof(ScalarType));
    ScalarType* zvalues = (ScalarType*) malloc(n*n*sizeof(ScalarType));
    ScalarType* yvalues = (ScalarType*) malloc(n*n*sizeof(ScalarType));

    /* On input, Z is equal to the symmetric matrix A.  On output, its columns
     * contain the n orthonormal eigenvectors of A.  Force the storage order to
     * be column-major, since it doesn't matter on input, and it makes things
     * easier for the output. */
    MyMatrix<ScalarType> Z(n, n, zvalues, ColumnMajor);
    A.copy(Z);

    lapack.SYEV(jobz, uplo, n, Z, Z.stride, w, work, lwork, &info);
    if(info != 0) {
        fprintf(stderr, "sqrt: SYEV returned %d\n", info);
    }
    
    /* Y = sqrt(W) . Z^T */
    MyMatrix<ScalarType> Y(n, n, yvalues, ColumnMajor);
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
        double s = sqrt(w[i]);  // sqrt of ith eigenvalue
        for(int j = 0; j < n; j++)
            Y(i,j) = s * Z(j,i);
    }

    /* B = Z . Y */
    multiply(Z, Y, B);

    /* Free temporary storage */
    free(w);
    free(work);
    free(zvalues);
    free(yvalues);
}

/* Invert the square matrix A, storing the result in B.  It is okay if B points
 * to the same memory as A.  */
template<class ScalarType>
void inverse(const MyMatrix<ScalarType>& A, MyMatrix<ScalarType>& B) {
#ifdef OPSEC_DEBUG
    assert(A.values != NULL && A.m == A.n);
    assert(B.values != NULL && B.m == A.m && B.n == A.n && B.storage == A.storage);
#endif

    Teuchos::LAPACK<int,ScalarType> lapack;
    int n = A.n;
    int blocksize = lapack.ILAENV(1, "DGETRI", "", n);
    int lwork = n*blocksize;    // optimal work size
    int info;

    /* Allocate workspace */
    int* ipiv = (int*) malloc(n*sizeof(int));
    ScalarType* work = (ScalarType*) malloc(lwork*sizeof(ScalarType));

    /* On input to GETRF, B is equal to the matrix A.  On output it contains
     * the factors L and U from the LU factorization A = P*L*U.  The pivot
     * indices determining the permutation matrix P are stored in ipiv. */
    if(B.values != A.values)
        A.copy(B);

    /* Perform LU factorization on the matrix A. */
    lapack.GETRF(n, n, B, B.stride, ipiv, &info);
    if(info != 0) {
        fprintf(stderr, "inverse: GETRF returned %d\n", info);
    }

    /* Use the LU factorization to compute the inverse of A. */
    lapack.GETRI(n, B, B.stride, ipiv, work, lwork, &info);
    if(info != 0) {
        fprintf(stderr, "inverse: GETRI returned %d\n", info);
    }

    /* Free workspace */
    free(ipiv);
    free(work);
}

#endif // MYMATRIX_HPP
