#ifndef EIG_H
#define EIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include "opsec.h"

/* Serial eigenvalue computation.  A is a real array of length n*n representing
 * an n-by-n symmetric matrix, evals is an array of length nev, and evecs is an
 * array of length n*nev.  ncv denotes the number of working vectors to use and
 * must lie between nev and n.  Returns the number of eigenvalues successfully
 * computed. */
int eig(int n, int nev, int ncv, real *A, real *evals, real *evecs);

/* Parallel eigenvalue computation.  A is a real array of length nloc*n
 * representing an nloc-by-n block of a symmetric matrix, stored in row-major
 * order.  evals is an array of length nev, and evecs is an array of length
 * nloc*nev.  ncv denotes the number of working vectors to use and must lie
 * between nev and n.  Returns the number of eigenvalues successfully computed. */
int peig(int n, int nloc, int nev, int ncv, real tol, real *A, real *evals, real *evecs);

#ifdef __cplusplus
}
#endif

#endif // EIG_H
