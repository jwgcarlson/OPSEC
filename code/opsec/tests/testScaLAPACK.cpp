/* testPdgemr2d
 *
 * Initializes a matrix A and a vector x on the root process, redistributes
 * them on a 2-D process grid, computes their product, then gathers the result
 * back to the root process. */

#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <cassert>
#include <cstdio>
#include <cstdlib>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "cscalapack.h"
#include "slp.h"

#if 0
/* BLACS context */
struct Context {
    int ictxt;      // BLACS context ID
    int nprow;      // number of rows in process grid
    int npcol;      // number of columns in process grid
    int myrow;      // my row
    int mycol;      // my column
};

Context create_context(int numrows, int numcols) {
    Context context;
    Cblacs_get(0, 0, &context.ictxt);
    Cblacs_gridinit(&context.ictxt, "Row", numrows, numcols);
    Cblacs_gridinfo(context.ictxt, &context.nprow, &context.npcol, &context.myrow, &context.mycol);
    if(context.myrow >= numrows || context.mycol >= numcols)
        context.myrow = context.mycol = -1;
    return context;
}

/* BLACS array descriptor */
struct ArrayDesc {
    int dtype;  // = 1 for dense matrices
    int ictxt;  // BLACS context handle
    int m;      // number of rows
    int n;      // number of columns
    int mb;     // row blocking factor
    int nb;     // column blocking factor
    int rsrc;   // = 0
    int csrc;   // = 0
    int lld;    // leading dimension of local array
};

static inline int mynumroc(int n, int nb, int iproc, int isrcproc, int nprocs) {
    return numroc_(&n, &nb, &iproc, &isrcproc, &nprocs);
}

static ArrayDesc create_descriptor(const Context& context,
                                   int m, int n, int mb, int nb,
                                   int rsrc = 0, int csrc = 0, int lld = 0)
{
    ArrayDesc desc;
    desc.dtype = 1;
    desc.ictxt = context.ictxt;
    desc.m = m;
    desc.n = n;
    desc.mb = mb;
    desc.nb = nb;
    desc.rsrc = rsrc;
    desc.csrc = csrc;
    desc.lld = std::max(lld, mynumroc(m, mb, context.myrow, rsrc, context.nprow));
    return desc;
}

//static inline int mydescinit(ArrayDesc* desc, int m, int n, int mb, int nb, int irsrc, int icsrc, const Context& context, int lld) {
//    int info;
//    printf("mydescinit: irsrc = %d, nprow = %d\n", irsrc, context.nprow);
//    descinit_((int*) desc, &m, &n, &mb, &nb, &irsrc, &icsrc, (int*) &context.ictxt, &lld, &info);
//    if(info != 0)
//        fprintf(stderr, "mydescinit: info = %d\n", info);
//    return info;
//}

static inline int myindxl2g(int iloc, int nb, int iproc, int isrcproc, int nprocs) {
    iloc += 1;
    return indxl2g_(&iloc, &nb, &iproc, &isrcproc, &nprocs) - 1;
}

static inline void mypdgemv(const char* transa, int m, int n,
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

extern "C" void Cpdgemr2d(int m, int n, double* a, int ia, int ja, int* desca, double* b, int ib, int jb, int* descb, int ictxt);
static inline void mypdgemr2d(int m, int n,
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

int main(int argc, char* argv[]) {
    /* Initialize MPI */
#ifdef OPSEC_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    /* Set up global BLACS context encompassing all processes */
    int nprocs, me;
    Cblacs_pinfo(&me, &nprocs);
    Context gcontext = create_context(nprocs, 1);
    me = gcontext.myrow;

    /* Decide on a good 2-D process grid for computation */
    int nprow = nprocs;
    int npcol = 1;
    while(nprow > npcol && (nprow % 2) == 0) {
        nprow /= 2;
        npcol *= 2;
    }

    /* Set up a new BLACS context for this 2-D process grid */
    Context pcontext = create_context(nprow, npcol);

    /* Finally, set up a BLACS context on the root process only for data I/O */
    Context rcontext = create_context(1, 1);
    if(me == 0)
        assert(rcontext.myrow == 0 && rcontext.mycol == 0);

    int m = 5;
    int n = 5;
    int mb = 2;
    int nb = 2;

    /* For the compute context, get the number of local rows and columns of the matrix A */
    int mloc = mynumroc(m, mb, pcontext.myrow, 0, pcontext.nprow);
    int nloc = mynumroc(n, nb, pcontext.mycol, 0, pcontext.npcol);

    printf("me = %d, myrow = %d, mycol = %d, mloc = %d, nloc = %d\n", me, pcontext.myrow, pcontext.mycol, mloc, nloc);

    /* Define array descriptors for local access of global arrays in the compute context */
    ArrayDesc descA = create_descriptor(pcontext, m, n, mb, nb, 0, 0, mloc);
    ArrayDesc descx = create_descriptor(pcontext, m, 1, mb, 1, 0, 0, mloc);
    ArrayDesc descy = create_descriptor(pcontext, m, 1, mb, 1, 0, 0, mloc);

    /* Allocate local storage */
    double* A = (double*) malloc(mloc*nloc*sizeof(double));
    double* x = (double*) malloc(mloc*sizeof(double));
    double* y = (double*) malloc(mloc*sizeof(double));

    /* Define array descriptor for a single column of A (only applicable on root process) */
    ArrayDesc desccolumn = create_descriptor(rcontext, m, 1, m, 1, 0, 0, m);

    /* Allocate memory for a full column on root process */
    double* column = NULL;
    if(me == 0)
        column = (double*) malloc(m*sizeof(double));

    /* For each column of A: set values on root process, then redistribute to compute context */
    for(int j = 0; j < n; j++) {
        if(me == 0) {
            /* Set values of column j of matrix A (or read them from file) */
            for(int i = 0; i < m; i++)
                column[i] = i - j;
        }

        /* Redistribute column to compute context */
        mypdgemr2d(m, 1, column, 0, 0, desccolumn, A, 0, j, descA, gcontext);
    }

    /* Set values of x on root process, then redstribute to compute context */
    if(me == 0) {
        for(int i = 0; i < m; i++)
            column[i] = (i % 2);
    }
    mypdgemr2d(m, 1, column, 0, 0, desccolumn, x, 0, 0, descx, gcontext);

    /* Multiply y = A*x */
    mypdgemv("N", m, n, 1.0, A, 0, 0, descA, x, 0, 0, descx, 1, 0.0, y, 0, 0, descy, 1);

    Cblacs_barrier(gcontext.ictxt, "A");

    /* Redistribute y back to root process, and print */
    mypdgemr2d(m, 1, y, 0, 0, descy, column, 0, 0, desccolumn, gcontext);
    if(me == 0) {
        for(int i = 0; i < m; i++)
            printf("y[%d] = %g\n", i, column[i]);
    }

    if(me == 0)
        free(column);
    free(A);
    free(x);
    free(y);
    Cblacs_gridexit(0);
    MPI_Finalize();
    return 0;
}
