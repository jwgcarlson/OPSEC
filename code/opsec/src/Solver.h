#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

#include "opsec.h"
#include "slp.h"

/* Abstract class for solving for the eigenvalues and eigenvectors of a
 * symmetric matrix. */
class Solver {
public:
    Solver();
    virtual ~Solver();

    virtual const slp::Context* GetContext() const = 0;

    /* Solve for eigenvalues and eigenvectors.  'nev' is the number requested.
     * Return the number actually obtained. */
    virtual int Solve(int nev) = 0;

    /* Retreive eigenvalues and eigenvectors after calling Solve(). */
    virtual real GetEigenvalue(int i) const = 0;
    virtual slp::Matrix<real> GetEigenvector(int i) const = 0;
};

/* Abstract class for computing selected values of a matrix. */
class MatrixFactory {
public:
    MatrixFactory() {}
    virtual ~MatrixFactory() {}

    /* Compute matrix elements A_{ij} for i in 'rows' and j in 'cols'.  Place
     * them in the column-major array 'values' with leading dimension 'lld'. */
    virtual void ComputeMatrixValues(
            const std::vector<int>& rows,
            const std::vector<int>& cols,
            real* values,
            int lld) = 0;
};

#endif // SOLVER_H
