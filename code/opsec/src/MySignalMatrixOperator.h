#ifndef MYSIGNALMATRIXOPERATOR_H
#define MYSIGNALMATRIXOPERATOR_H

#include <algorithm>
#include <vector>

#include <AnasaziOperator.hpp>
#include <Teuchos_BLAS.hpp>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#  include <Teuchos_RawMPITraits.hpp>
#endif

#include "Cell.h"
#include "Model.h"
#include "MyMultiVec.h"
#include "XiFunc.h"
#include "sig.h"


template<class ScalarType>
class MySignalMatrixOperator : public Anasazi::Operator<ScalarType> {
public:
    typedef MyMultiVec<ScalarType> MV;

    /* Constructor */
    MySignalMatrixOperator(int coordsys_, int N1_, int N2_, int N3_,
            int Ncells_, const Cell* cells_, XiFunc xi_, Survey* survey)
    {
        coordsys = coordsys_;
        N1 = N1_;
        N2 = N2_;
        N3 = N3_;

        Ncells = Ncells_;
        cells = cells_;

        xi = xi_;

        nprocs = 1;
        me = 0;
#ifdef OPSEC_USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
        n.resize(nprocs);
        d.resize(nprocs);
        for(int p = 0; p < nprocs; p++) {
            n[p] = (Ncells/nprocs) + (p < (Ncells % nprocs));
            d[p] = p*(Ncells/nprocs) + std::min(p, Ncells % nprocs);
        }

        opsec_info("Computing signal matrix components...\n");

        std::vector<int> rows(n[me]), cols(Ncells);
        for(int aloc = 0; aloc < n[me]; aloc++)
            rows[aloc] = aloc + d[me];
        for(int b = 0; b < Ncells; b++)
            cols[b] = b;

        S = (real*) opsec_malloc(n[me]*Ncells*sizeof(real));
        ComputeSignalMatrixBlock(Ncells, rows, cols, S, n[me],
                                 xi, survey, coordsys, cells, N1, N2, N3);
    }

    /* Destructor */
    ~MySignalMatrixOperator() {
        free(S);
    }

    void Apply(const Anasazi::MultiVec<ScalarType>& x, Anasazi::MultiVec<ScalarType>& y) const
    {
        assert(x.GetVecLength() == Ncells && y.GetVecLength() == Ncells);
        assert(x.GetNumberVecs() == y.GetNumberVecs());

        /* Cast to my types */
        const MyMultiVec<ScalarType>* px = dynamic_cast<const MyMultiVec<ScalarType>*>(&x);
        MyMultiVec<ScalarType>* py = dynamic_cast<MyMultiVec<ScalarType>*>(&y);
        assert(px != NULL && py != NULL);

        std::vector<ScalarType> xfull(Ncells);  // one full vector
        int numvecs = x.GetNumberVecs();
        for(int j = 0; j < numvecs; j++) {
            /* Copy full vector j of x into xfull */
            MPI_Allgatherv((void*) px->vector(j),
                           n[me],
                           Teuchos::RawMPITraits<ScalarType>::type(),
                           &xfull[0],
                           (int*) &n[0],
                           (int*) &d[0],
                           Teuchos::RawMPITraits<ScalarType>::type(),
                           MPI_COMM_WORLD);

            /* Apply local signal matrix to xfull */
            blas.GEMV(Teuchos::NO_TRANS,// trans
                      n[me],            // m
                      Ncells,           // n
                      1,                // alpha
                      S,                // a
                      Ncells,           // lda
                      &xfull[0],        // x
                      1,                // incx
                      0,                // beta
                      py->vector(j),    // y
                      1);               // incy
        }
    }

private:
    /* Coordinate system: either CoordSysCartesian or CoordSysSpherical */
    int coordsys;

    /* Number of cell divisions in each coordinate direction */
//    union { double N1; double Nx; double Nr; };
//    union { double N2; double Ny; double Nmu; };
//    union { double N3; double Nz; double Nphi; };
    int N1, N2, N3;

    int Ncells;         // a.k.a. N, total matrix/vector size
    const Cell* cells;

    XiFunc xi;

    Teuchos::BLAS<int, ScalarType> blas;        // BLAS interface

//    std::vector<ScalarType> Q;  // cache of pre-computed values

    int nprocs;         // number of processes
    int me;             // my process rank

    std::vector<int> n; // local matrix/vector size for each process:
                        //   n[p] = (N/nprocs) + (p < (N % nprocs))
    std::vector<int> d; // first element assigned to each process:
                        //   d[p] = p*(N/nprocs) + std::min(p, N % nprocs)
                        //        = n[0] + n[1] + ... + n[p-1]

    ScalarType* S;      // my local portion of the signal matrix
};

#endif // MYSIGNALMATRIXOPERATOR_H
