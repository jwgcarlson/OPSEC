#include <algorithm>
#include <cassert>
#include <cstring>
#include <ctime>

#include <mpi.h>

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBasicSort.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "AnasaziSolver.h"
#include "MyAnasazi.h"
#include "opsec.h"
#include "slp.h"

#if 0
/* Return n divided by p, rounded up to nearest integer */
static inline int divup(int n, int p) {
    return (n + p - 1)/p;
}
#endif

#if 0
void old() {
    /* Number of eigenvalues to compute */
    int nev = Nmodes;

    /* Block size to use for block Krylov-Schur method */
    int block_size = 1;

    /* Number of blocks to use */
    int num_blocks = std::min(2*nev/block_size, Ncells);

    /* Convergence tolerance */
    double tol = 1e-5;

    /* Maximum number of restarts */
    int max_restarts = 1000;

    /* Which eigenvalues to solve for (LM = largest magnitude) */
    std::string which = "LM";

    /* Verbosity level */
    int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary
                  + Anasazi::TimingDetails + Anasazi::StatusTestDetails;

    /* Construct parameter list for Eigensolver */
    Teuchos::ParameterList params;
    params.set("Verbosity", verbosity);
    params.set("Block Size", block_size);
    params.set("Num Blocks", num_blocks);
    params.set("Maximum Restarts", max_restarts);
    params.set("Convergence Tolerance", tol);
    params.set("Which", which);

    typedef real ST;
    typedef Anasazi::MultiVec<ST> MV;
    typedef Anasazi::Operator<ST> OP;
    typedef Anasazi::MultiVecTraits<ST,MV> MVT;
    typedef Anasazi::OperatorTraits<ST,MV,OP> OPT;

    /* Construct initial vector */
    Teuchos::RCP<MyMultiVec<real> > ivec = Teuchos::rcp(new MyMultiVec<real>(Ncells, block_size));
    ivec->MvRandom();

    if(me == 0)
        opsec_info("Defining signal matrix operator...\n");

    /* Construct operator */
    Teuchos::RCP<MySignalMatrixOperator<real> > op
        = Teuchos::rcp(new MySignalMatrixOperator<real>(coordsys, N1, N2, N3, Ncells, cells, xi, survey, epsrel, epsabs));

    if(me == 0)
        opsec_info("Initializing eigenproblem...\n");

    /* Create the eigenproblem */
    Teuchos::RCP<Anasazi::BasicEigenproblem<ST, MV, OP> > problem =
        Teuchos::rcp(new Anasazi::BasicEigenproblem<ST, MV, OP>(op, ivec));
    problem->setHermitian(true);
    problem->setNEV(nev);
    bool problem_okay = problem->setProblem();
    if(!problem_okay && me == 0) {
        opsec_error("klt: error setting eigenproblem\n");
        opsec_abort(1);
    }

    if(me == 0)
        opsec_info("Initializing eigensolver...\n");

    /* Initialize the block-Arnoldi solver */
    Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> solmgr(problem, params);

    if(me == 0)
        opsec_info("Solving eigenproblem...\n");

    /* Solve the problem to the specified tolerances or length */
    Anasazi::ReturnType retcode = solmgr.solve();
    if(retcode != Anasazi::Converged && me == 0) {
        opsec_error("klt: SolverManager failed to converge\n");
    }

    if(me == 0)
        opsec_info("Reading solution...\n");

    /* Get the eigenvalues and eigenvectors from the eigenproblem */
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    std::vector<Anasazi::Value<ST> > evals = sol.Evals;
    Teuchos::RCP<MV> mvevecs = sol.Evecs;
    MyMultiVec<ST>* evecs = dynamic_cast<MyMultiVec<ST>*>(mvevecs.get());
    std::vector<int> index = sol.index;
    int nconv = sol.numVecs;
}
#endif

AnasaziSolver::AnasaziSolver(int n, MatrixFactory& matfact)
    : n(n), nconv(0)
{
    /* Initialize BLACS process grid */
    int nprocs, me;
    Cblacs_pinfo(&me, &nprocs);
    pcontext = new slp::Context(nprocs, 1);

    nb = divup(n, nprocs);

    /* Initialize column vector descriptor */
    xdesc = pcontext->new_descriptor(n, 1, nb, 1);

    /* Instantiate operator */
    op = Teuchos::rcp(new MyOperator(pcontext, n, matfact, nb));

    /* Set default Anasazi parameters */
    block_size = 1;
    num_blocks = 0;
    tol = 1e-5;
    max_restarts = 1000;
    which = "LM";               // LM = largest magnitude
    verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary
                  + Anasazi::TimingDetails + Anasazi::StatusTestDetails;
}

AnasaziSolver::~AnasaziSolver() {
    delete xdesc;

    pcontext->exit();
    delete pcontext;
}

int AnasaziSolver::Solve(int nev) {
    int me, nprocs;
    Cblacs_pinfo(&me, &nprocs);

    /* Construct initial vector */
    Teuchos::RCP<MyMultiVec> ivec = Teuchos::rcp(new MyMultiVec(pcontext, n, block_size, nb));
    ivec->MvRandom();

    if(me == 0)
        opsec_info("Initializing eigenproblem...\n");

    /* Create the eigenproblem */
    Teuchos::RCP<Anasazi::BasicEigenproblem<ST, MV, OP> > problem =
        Teuchos::rcp(new Anasazi::BasicEigenproblem<ST, MV, OP>(op, ivec));
    problem->setHermitian(true);
    problem->setNEV(nev);
    bool problem_okay = problem->setProblem();
    if(!problem_okay && me == 0) {
        opsec_error("AnasaziSolver::Solve: error setting eigenproblem\n");
        opsec_abort(1);
    }

    if(me == 0)
        opsec_info("Initializing eigensolver...\n");

    /* Choose reasonable number of blocks for computation */
    if(num_blocks <= 0)
        num_blocks = divup(2*nev, block_size);
    num_blocks = std::min(n, num_blocks);

    /* Construct parameter list for Eigensolver */
    Teuchos::ParameterList params;
    params.set("Block Size", block_size);
    params.set("Num Blocks", num_blocks);
    params.set("Convergence Tolerance", tol);
    params.set("Maximum Restarts", max_restarts);
    params.set("Which", which);
    params.set("Verbosity", verbosity);

    /* Initialize the block-Arnoldi solver */
    Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> solmgr(problem, params);

    if(me == 0)
        opsec_info("Solving eigenproblem...\n");

    /* Solve the problem to the specified tolerances or length */
    Anasazi::ReturnType retcode = solmgr.solve();
    if(retcode != Anasazi::Converged && me == 0) {
        opsec_error("AnasaziSolver::Solve: SolverManager failed to converge\n");
    }

    if(me == 0)
        opsec_info("Reading solution...\n");

    /* Get the eigenvalues and eigenvectors from the eigenproblem */
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    std::vector<Anasazi::Value<ST> > complex_evals = sol.Evals;
    evecs = sol.Evecs;
    index = sol.index;
    nconv = sol.numVecs;
    for(int i = 0; i < nconv; i++) {
        assert(fabs(complex_evals[i].imagpart) < 1e-12);
        evals[i] = complex_evals[i].realpart;
    }
}

const slp::Context* AnasaziSolver::GetContext() const {
    return pcontext;
}

real AnasaziSolver::GetEigenvalue(int i) const {
    assert(0 <= i && i < nconv);
    return evals[i];
}

slp::Matrix<real> AnasaziSolver::GetEigenvector(int i) const {
    assert(0 <= i && i < nconv);
    assert(index[i] == 0);
    MyMultiVec* v = dynamic_cast<MyMultiVec*>(evecs.get());
    return slp::Matrix<real>(xdesc, v->vector(i));
}
