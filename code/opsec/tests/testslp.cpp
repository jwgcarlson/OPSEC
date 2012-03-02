
#include <cstdio>
#include <cstdlib>

#include <mpi.h>

#include "slp.h"

using namespace slp;

int main(int argc, char* argv[]) {
    /* Matrix dimensions */
    const int N = 50;
    const int NB = 16;

    int nprocs = 1, me = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    /* Create a roughly square process grid for computation */
    int nprow = nprocs;
    int npcol = 1;
    while(npcol < nprow && (nprow % 2) == 0) {
        nprow /= 2;
        npcol *= 2;
    }
    Context pcontext(nprow, npcol);

    /* Create a N-by-N matrix distributed in NB-by-NB blocks */
    Descriptor* Adesc = pcontext.new_descriptor(N, N, NB, NB);
    double* Avalues = (double*) malloc(Adesc->local_size() * sizeof(double));
    Matrix<double> A(Adesc, Avalues);

    /* Initialize A to a "shift by one" matrix */
    for(LocalMatrixIndex k(A); !k.done(); ++k) {
        int i = k.i, j = k.j;
        A[k] = (i+1 == j);
    }

    /* Create N-by-4 multi-vectors x and y */
    Descriptor* xdesc = pcontext.new_descriptor(N, 4, NB, 4);
    double* xvalues = (double*) malloc(xdesc->local_size() * sizeof(double));
    double* yvalues = (double*) malloc(xdesc->local_size() * sizeof(double));
    Matrix<double> x(xdesc, xvalues);
    Matrix<double> y(xdesc, yvalues);

    /* Initialize columns of x */
    for(LocalMatrixIndex k(x); !k.done(); ++k) {
        int i = k.i, j = k.j;
        x[k] = i + j;
    }

    /* Compute y = A*x */
    multiply(A, x, y);

    /* Create N-by-1 vector on root process only */
    Descriptor* zdesc = pcontext.new_descriptor(N, 1, N, 1);
    double* zvalues = (double*) malloc(zdesc->local_size() * sizeof(double));
    Matrix<double> z(zdesc, zvalues);

    /* Gather columns of y to root process one at a time, and print */
    for(int j = 0; j < 4; j++) {
        redistribute(N, 1, y, 0, j, z, 0, 0);
        if(me == 0) {
            for(int i = 0; i < N; i++)
                printf("y[%d] = %g\n", i, z[i]);
        }
    }

    /* Clean up */
    delete Adesc;
    delete xdesc;
    delete zdesc;
    free(Avalues);
    free(xvalues);
    free(yvalues);
    free(zvalues);
//    pcontext.exit();
    Cblacs_gridexit(0);
    MPI_Finalize();
}
