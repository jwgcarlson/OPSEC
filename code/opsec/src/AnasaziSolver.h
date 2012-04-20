#ifndef ANASAZISOLVER_H
#define ANASAZISOLVER_H

#include <string>

#include <AnasaziMultiVec.hpp>
#include <AnasaziOperator.hpp>
#include <Teuchos_RCP.hpp>

#include "Solver.h"

class AnasaziSolver : public Solver {
public:
    AnasaziSolver(int n, MatrixFactory& matfact);
    ~AnasaziSolver();

    const slp::Context* GetContext() const;

    /* Solve for eigenvalues and eigenvectors.  nev is the number requested.
     * Return the number actually obtained. */
    int Solve(int nev);

    /* Retreive eigenvalues and eigenvectors after calling Solve(). */
    real GetEigenvalue(int i) const;
    slp::Matrix<real> GetEigenvector(int i) const;

private:
    typedef real ST;
    typedef Anasazi::MultiVec<ST> MV;
    typedef Anasazi::Operator<ST> OP;
    typedef Anasazi::MultiVecTraits<ST,MV> MVT;
    typedef Anasazi::OperatorTraits<ST,MV,OP> OPT;

    /* Global matrix size */
    int n;

    /* Row-wise blocking factor */
    int nb;

    /* Number of converged eigenvalues, after calling Solve(). */
    int nconv;

    /* Anasazi parameters */
    int block_size;     // block size for block Krylov-Schur method
    int num_blocks;     // number of blocks to use
    double tol;         // convergence tolerance
    int max_restarts;   // maximum number of restarts
    std::string which;  // which eigenvalue to choose (LM = largest magnitude)
    int verbosity;      // verbosity level

    /* Anasazi products */
    Teuchos::RCP<OP> op;
    Teuchos::RCP<MV> evecs;
    std::vector<ST> evals;
    std::vector<int> index;

    /* BLACS process grid */
    slp::Context* pcontext;

    /* Descriptor for a single eigenvector */
    slp::Descriptor* xdesc;
};

#endif // ANASAZISOLVER_H
