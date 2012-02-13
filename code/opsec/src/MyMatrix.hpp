#ifndef MYMATRIX_HPP
#define MYMATRIX_HPP

#include <cassert>

#include <Teuchos_BLAS.hpp>

/* Simple classes for treating regular arrays of numbers as vectors or matrices. */


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

    MyVector(int n_, int inc_, ScalarType* values_) {
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

    /* Copy this vector into a new vector.  The array newvalues must be able
     * to hold at least n elements, and must not overlap with this->values. */
    MyVector copy(ScalarType* newvalues) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && newvalues != NULL);
#endif
        COPY(n, values, inc, newvalues, 1);
        return MyVector(n, 1, newvalues);
    }


    MagnitudeType norm() const {
#ifdef OPSEC_DEBUG
        assert(values != NULL && n > 0);
#endif
        return NRM2(n, values, inc);
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
struct MyMatrix {

    /***** Member variables **************************************************/

    int m;              // number of rows
    int n;              // number of columns
    int storage;        // storage order, either ColumnMajor or RowMajor
    int stride;         // column or row stride, depending on storage order
    ScalarType* values;


    /***** Constructors *******************************************************/

    MyMatrix() {
        m = n = stride = 0;
        storage = -1;
        values = NULL;
    }

    /* Set stride to 0 to indicate packed storage, i.e.
     *   stride = m if storage == ColumnMajor,
     *   stride = n if storage == RowMajor. */
    MyMatrix(int m_, int n_, int storage_, int stride_, ScalarType* values_) {
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


    /* Create a view into a m1-by-n1 sub-matrix of this matrix. */
    MyMatrix block(int m1, int n1, int i1, int j1) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && m1 > 0 && n1 > 0 && 0 <= i1 && i1 + m1 <= m && 0 <= j1 && j1 + n1 <= n);
#endif
        return MyMatrix(m1, n1, storage, stride, &element(i1,j1));
    }

    /* Copy values into a new matrix.  The array newvalues must be able to hold
     * at least m*n elements, and must not overlap with this->values.  The new
     * matrix has packed storage, i.e. there is no padding after each row or
     * column. */
    MyMatrix copy(ScalarType* newvalues) {
#ifdef OPSEC_DEBUG
        assert(values != NULL && newvalues != NULL);
#endif
        if((storage == ColumnMajor && stride == m) || (storage == RowMajor && stride == n)) {
            /* Packed matrix storage, memcpy the entire array */
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
        return MyMatrix(m, n, storage, newstride, newvalues);
    }

    /* Return matrix transpose. */
    MyMatrix transpose() {
        int newstorage = (storage == ColumnMajor) ? RowMajor : ColumnMajor;
        return MyMatrix(m, n, newstorage, stride, values);
    }

    /***** Accessors **********************************************************/

    operator ScalarType*() {
        return values;
    }

    operator const ScalarType*() const {
        return values;
    }

    ScalarType& element(int i, int j) {
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

    MyVector<ScalarType> column(int j) {
        if(storage == ColumnMajor)
            return MyVector<ScalarType>(m, 1, &values[j*stride]);
        else
            return MyVector<ScalarType>(m, stride, &values[j]);
    }

    MyVector<ScalarType> row(int i) {
        if(storage == ColumnMajor)
            return MyVector<ScalarType>(n, stride, &values[i]);
        else
            return MyVector<ScalarType>(n, 1, &values[i*stride]);
    }

};


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

/* Matrix-matrix product */
template<class ScalarType>
void multiply(const MyMatrix<ScalarType>& A, const MyMatrix<ScalarType>& B, MyMatrix<ScalarType>& C) {
#ifdef OPSEC_DEBUG
    assert(A.values != NULL && B.values != NULL && C.values != NULL && C.m == A.m && A.n == B.m && B.n == C.n);
#endif
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::ETransp transa = (A.storage == ColumnMajor) ? Teuchos::NO_TRANS : Teuchos::TRANS;
    Teuchos::ETransp transb = (A.storage == ColumnMajor) ? Teuchos::NO_TRANS : Teuchos::TRANS;
    blas.GEMM(transa, transb, C.m, C.n, A.n, 1, A, A.stride, B. B.stride, 0, C, C.stride);
}

#endif // MYMATRIX_HPP
