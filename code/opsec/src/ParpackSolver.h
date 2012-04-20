#ifndef PARPACKSOLVER_H
#define PARPACKSOLVER_H

#include "Solver.h"

class ParpackSolver : public Solver {
public:
    ParpackSolver(int n, MatrixFactory& matfact);
    ~ParpackSolver();

    const slp::Context* GetContext() const;
    slp::Matrix<real> GetMatrix() const;

    /* Solve for eigenvalues and eigenvectors.  nev is the number requested.
     * Return the number actually obtained. */
    int Solve(int nev);

    real GetEigenvalue(int i) const;
    slp::Matrix<real> GetEigenvector(int i) const;

//    const std::vector<real>& GetEigenvalues() const;
//    slp::Matrix<real> GetEigenvectors() const;

private:
    /* Global matrix size */
    int n;

    /* BLACS process grid */
    slp::Context* pcontext;

    /* Symmetric matrix $A$, shared over all processes.  This matrix is
     * initialized in the constructor. */
    slp::Descriptor* Adesc;
    real* Avalues;

    /* Number of converged eigenvalues after call to Solve(). */
    int nconv;

    /* Array of eigenvalues after calling Solve().  Each process stores the
     * full array. */
    std::vector<real> evals;

    /* Descriptor for a length-n column vector shared over all processes. */
    slp::Descriptor* xdesc;

    /* Descriptor for a n-by-ncv matrix $B$ shared over all processes, used to
     * store working Ritz vectors and, after calling Solve(), the full set of
     * eigenvectors. */
    slp::Descriptor* Bdesc;

    /* Matrix of eigenvectors after calling Solve().  Data layout is specified
     * by Bdesc. */
    real* Bvalues;

    /* PARPACK parameters */
    int nloc;           // number of local rows for my process
    double tol;         // convergence tolerance
    int maxitr;         // maximum number of Arnoldi iterations
    int ncv;            // number of working Ritz vectors to use
};

#endif // PARPACKSOLVER_H
