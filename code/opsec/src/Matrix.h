#ifndef MATRIX_H
#define MATRIX_H

/* TODO:
 * - refactor operator+/-(), transpose(), etc. to take advantage of RVO */

/* A simple class for arbitrarily-sized matrices. */
template<typename T>
class Matrix {
public:
    /* Construct an empty matrix. */
    Matrix();

    /* Construct an (uninitialized) MxN matrix. */
    Matrix(int M, int N);

    /* Construct MxN matrix from C array, either by copying the elements or
     * taking direct ownership of the (malloc-allocated) memory */
    Matrix(int M, int N, T* b, bool copy = true);

    /* Copy constructor. */
    Matrix(const Matrix<T>& B);

    /* Assignment. */
    Matrix& operator=(const Matrix<T>& B);

    /* Destructor. */
    ~Matrix();

    /* Set all elements to zero. */
    void zero();

    /* Resize the matrix.  Previously stored values are lost. */
    void resize(int M, int N);

    /* Matrix element accessors. */
    const T&  operator()(int i, int j = 0) const;
    T& operator()(int i, int j = 0);

    /* Unchecked data element accessors. */
    const T&  operator[](int i) const;
    T& operator[](int i);

    /* Auto-cast to scalar for 1-by-1 matrices */
    operator T() const;

public:         // internal data is explicitly public
    int M, N;
    T* data;
};

template<typename T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> dot(const Matrix<T>& A, const Matrix<T>& B);

template<typename T>
Matrix<T> transpose(const Matrix<T>& A);

template<typename T>
Matrix<T> identity(int N);


/***** Inline implementation *****/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "opsec.h"

template<typename T>
inline Matrix<T>::Matrix() {
    M = N = 0;
    data = NULL;
}

template<typename T>
inline Matrix<T>::Matrix(int m, int n) {
    M = m;
    N = n;
    data = (T*)malloc(M*N*sizeof(T));
}

template<typename T>
inline Matrix<T>::Matrix(int m, int n, T* b, bool copy) {
    M = m;
    N = n;
    if(b == NULL)
        data = (T*)malloc(M*N*sizeof(T));
    else if(!copy)
        data = b;
    else {
        data = (T*)malloc(M*N*sizeof(T));
        memcpy(data, b, M*N*sizeof(T));
    }
}

template<typename T>
inline Matrix<T>::Matrix(const Matrix<T>& B) {
    M = B.M;
    N = B.N;
    data = (T*)malloc(M*N*sizeof(T));
    if(B.data)
        memcpy(data, B.data, M*N*sizeof(T));
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& B) {
    free(data);
    M = B.M;
    N = B.N;
    data = (T*)malloc(M*N*sizeof(T));
    if(B.data)
        memcpy(data, B.data, M*N*sizeof(T));
    return *this;
}

template<typename T>
inline Matrix<T>::~Matrix() {
    free(data);
}

template<typename T>
inline void Matrix<T>::zero() {
    memset(data, 0, M*N*sizeof(T));
}

template<typename T>
inline void Matrix<T>::resize(int m, int n) {
    M = m;
    N = n;
    data = (T*)realloc(data, M*N*sizeof(T));
}

template<typename T>
inline const T& Matrix<T>::operator()(int i, int j) const {
    assert(0 <= i && i < M && 0 <= j && j < N);
    return data[N*i + j];
}

template<typename T>
inline T& Matrix<T>::operator()(int i, int j) {
    assert(0 <= i && i < M && 0 <= j && j < N);
    return data[N*i + j];
}

template<typename T>
inline const T& Matrix<T>::operator[](int i) const {
    return data[i];
}

template<typename T>
inline T& Matrix<T>::operator[](int i) {
    return data[i];
}

template<typename T>
inline Matrix<T>::operator T() const {
    assert(M == 1 && N == 1 && data != NULL);
    return data[0];
}

/* Fortran expects multi-dimensional arrays to be stored in column-major order,
 * but in C we choose to store them in row-major order.  Thus, in C, A is an
 * m-by-k matrix and B is a k-by-n matrix, but when passed to a Fortran routine,
 * they appear as k-by-m and n-by-k matrices, respectively.  So to obtain A*B
 * we actually compute transpose(transpose(B)*transpose(A)), where all three
 * transposes are taken care of naturally by the Fortran-to-C interchange. */
template<typename T>
inline Matrix<T> dot(const Matrix<T>& A , const Matrix<T>& B) {
    int m = A.M;
    int k = A.N;
    assert(B.M == k);
    int n = B.N;
    real* a = const_cast<real*>(&A(0,0));
    real* b = const_cast<real*>(&B(0,0));
    real* c = (real*)malloc(m*n*sizeof(real));
    char notrans = 'N';
    real alpha = 1.0, beta = 0.0;
    blas_gemm(&notrans, &notrans, &n, &m, &k, &alpha, b, &n, a, &k, &beta, c, &n);
    return Matrix<real>(m, n, c, false);
}

template<typename T>
inline Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
    int M = A.M, N = A.N;
    assert(B.M == M && B.N == N);
    Matrix<T> C(M, N);
    int i, j;
    #pragma omp parallel for private(i,j)
    for(i = 0; i < M; i++)
        for(j = 0; j < N; j++)
            C(i,j) = A(i,j) + B(i,j);
    return C;
}

template<typename T>
inline Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
    int M = A.M, N = A.N;
    assert(B.M == M && B.N == N);
    Matrix<T> C(M, N);
    int i, j;
    #pragma omp parallel for private(i,j)
    for(i = 0; i < M; i++)
        for(j = 0; j < N; j++)
            C(i,j) = A(i,j) - B(i,j);
    return C;
}

template<typename T>
inline Matrix<T> transpose(const Matrix<T>& A) {
    int M = A.M, N = A.N;
    Matrix<T> At(N, M);
    int i, j;
    #pragma omp parallel for private(i,j)
    for(i = 0; i < M; i++)
        for(j = 0; j < N; j++)
            At(j,i) = A(i,j);
    return At;
}

template<typename T>
inline Matrix<T> identity(int N) {
    Matrix<T> I(N, N);
    int i, j;
//    #pragma omp parallel for private(i,j)
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            I(i,j) = (i == j);
    return I;
}

#endif // MATRIX_H
