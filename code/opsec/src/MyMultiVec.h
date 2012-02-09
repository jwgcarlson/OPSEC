#ifndef MYMULTIVEC_H
#define MYMULTIVEC_H

/* Concrete implementation of Anasazi::MultiVec<ScalarType> class.
 * - Employs distributed data model using MPI, with each vector split across
 *   processes in the natural way.  For a vector of length N distributed over
 *   k processes, the process with rank p is assigned
 *       Nloc = (N/k) + (p < (N % k))
 *   rows of the vector.  That is, if k evenly divides N, then each process
 *   gets N/k rows; otherwise, if the remainder is r = (N % k), then the first
 *   r processes get floor(N/k)+1 rows, while the remaining processes get
 *   floor(N/k) rows each.  [This is a "natural" data distribution because,
 *   amongst other reasons, it makes it very easy to decide which process a
 *   given row i of the vector belongs to: p = floor(i*k/N).]
 * - For each process, local data is further broken down into "data blocks",
 *   which is an abstraction for a contiguous block of data in memory.  The
 *   blocks need not all be the same size.  This additional complication is
 *   introduced to help optimize the most common use cases for MultiVec's.
 * - Overall the data distribution looks schematically as follows:
 *                                 Vectors
 *                  0  1  2  3  4  5  6  7  8  9  10 11 12 13
 *                 |  :  :  :  |  :  :  :  :  |  :  :  :  :  |
 *                 +-----------+--------------+--------------+
 *                 |  :  :  :  |  :  :  :  :  |  :  :  :  :  |
 *     Process 0   |  block 0  |   block 1 :  |  : block 2   |
 *                 |  :  :  :  |  :  :  :  :  |  :  :  :  :  |
 *                 +-----------+--------------+--------------+
 *                 |  :  :  :  |  :  :  :  :  |  :  :  :  :  |
 *     Process 1   |  block 0  |   block 1 :  |  : block 2   |
 *                 |  :  :  :  |  :  :  :  :  |  :  :  :  :  |
 *                 +-----------+--------------+--------------+
 *       ...       |           |              |              |
 *
 *   Each block is a contiguous region of local memory, representing m vectors
 *   each of local length n.
 */

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <AnasaziMultiVec.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#  include <Teuchos_RawMPITraits.hpp>
#endif

#include "abn.h"

template<typename T>
static inline int sgn(T x) {
    if(x < 0)      return -1;
    else if(x > 0) return +1;
    else           return 0;
}

/* Decide if two arrays overlap in memory.  Compare the sign of (&b[0]-&a[0])
 * and (&b[0]-&a[na]); if they differ, then memory region a straddles the start
 * of memory region b.  And vice versa to test if b straddles a. */
template<typename T>
static inline bool array_overlap(const T* a, int na, const T* b, int nb) {
    int x1 = sgn(b - a) * sgn(b - a - na);
    int x2 = sgn(a - b) * sgn(a - b - nb);
    return (x1 == -1) || (x2 == -1);
}


/* MyMultiVec<ScalarType>
 *
 * Distributed-memory multi-vector class, using MPI. */
template<class ScalarType>
class MyMultiVec : public Anasazi::MultiVec<ScalarType> {
public:

    /* Constructor */
    MyMultiVec(int length_, int numvecs_ = 0) {
        length = length_;
        numvecs = numvecs_;
        assert(length > 0);

        nprocs = 1;
        me = 0;
#ifdef OPSEC_USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

        mylength = (length/nprocs) + (me < (length % nprocs));
        mystart = me*(length/nprocs) + std::min(me, length % nprocs);

        if(numvecs > 0) {
            /* Allocate local memory, and assign all vectors to a single data block */
            DataBlock b(mylength, numvecs);
            blocks.push_back(b);
            block_widths.push_back(numvecs);
            block_starts.push_back(&b[0]);
            for(int i = 0; i < numvecs; i++) {
                VectorInfo info = { 0, i };
                index_to_block_map.push_back(info);
            }
        }
    }

    /* Destructor */
    ~MyMultiVec() {
        /* Cleanup should be automatic */
    }


    /***** Input/output *******************************************************/

    /* Read multi-vector from a .abn file.  All sizes must match.  The binary
     * data is arranged in column-major format, i.e. one full vector after
     * another, so it must be transposed to match the distributed storage
     * layout of the multi-vector. */
    void ReadFromFile(const char* filename) {
        int n;  // length of vector
        int m;  // number of vectors
        real* modes = NULL;
        int nread;

        int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

        if(me == 0) {
            /* Read in spectrum from file */
            FILE* fp = fopen(modefile, "r");
            if(fp == NULL) {
                fprintf(stderr, "ReadModesHeader: could not open file '%s'\n", modefile);
                return NULL;
            }

            size_t n, size;
            Config opts = cfg_new();
            int err = abn_read_header(fp, &n, &size, NULL, NULL, opts);
            int Nmodes = cfg_get_int(opts, "Nmodes");
            int Ncells = cfg_get_int(opts, "Ncells");
            cfg_destroy(opts);

            if(err || Nmodes <= 0 || Ncells <= 0 || n != Nmodes*Ncells || size != sizeof(double)) {
                fprintf(stderr, "ReadModesHeader: error reading '%s'\n", modefile);
                fclose(fp);
                return NULL;
            }

            if(fp != NULL) {
                modes = (real*) opsec_malloc(Nmodes*Ncells*sizeof(real));
                if(modes == NULL)
                    fprintf(stderr, "ReadModes: could not allocate %zd bytes of memory\n", Nmodes*Ncells*sizeof(real));
                else {
                    if(fread(modes, sizeof(real), Nmodes*Ncells, fp) != Nmodes*Ncells) {
                        fprintf(stderr, "ReadModes: error reading data from '%s'\n", modefile);
                        free(modes);
                        modes = NULL;
                    }
                }
                fclose(fp);
            }
        }

#ifdef OPSEC_USE_MPI
        /* Broadcast spectrum to other processes */
        if(usempi) {
            MPI_Bcast(&Nmodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&Ncells, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(Nmodes > 0) {
                if(me != 0)
                    modes = (real*) opsec_malloc(Nmodes*Ncells*sizeof(real));
                MPI_Bcast(modes, Nmodes*Ncells, REAL_MPI_TYPE, 0, MPI_COMM_WORLD);
            }
        }
#endif

        if(Nmodes_) *Nmodes_ = Nmodes;
        if(Ncells_) *Ncells_ = Ncells;
        return modes;
    }


    /***** Accessors **********************************************************/

    /* Access elements by local row index i */
    ScalarType& local(int i, int j) {
        assert(0 <= i && i < mylength && 0 <= j && j < numvecs);
        /* Find block that owns vector number j */
        int b = index_to_block_map[j].b;
        int k = index_to_block_map[j].k;
        return block_starts[b][k*mylength + i];
    }

    const ScalarType& local(int i, int j) const {
        assert(0 <= i && i < mylength && 0 <= j && j < numvecs);
        int b = index_to_block_map[j].b;
        int k = index_to_block_map[j].k;
        return block_starts[b][k*mylength + i];
    }

    /* Access elements by global row index i (read-only) */
    ScalarType global(int i, int j) {
        assert(0 <= i && i < length && 0 <= j && j < numvecs);
        ScalarType x;
        int owner = i*nprocs/length;    // process that owns row i of the multi-vector
        if(me == owner)
            x = local(i - mystart, j);
#ifdef OPSEC_USE_MPI
        MPI_Bcast(&x, 1, Teuchos::RawMPITraits<ScalarType>::type(), owner, MPI_COMM_WORLD);
#endif
        return x;
    }

    /* Access local vectors */
    ScalarType* vector(int j) {
        assert(0 <= j && j < numvecs);
        /* Find block that owns vector number j */
        int b = index_to_block_map[j].b;
        int k = index_to_block_map[j].k;
        return block_starts[b] + k*mylength;
    }

    const ScalarType* vector(int j) const {
        assert(0 <= j && j < numvecs);
        /* Find block that owns vector number j */
        int b = index_to_block_map[j].b;
        int k = index_to_block_map[j].k;
        return block_starts[b] + k*mylength;
    }


    /***** Anasazi::MultiVec<ScalarType> interface ****************************/

    /* Creation methods */

    MyMultiVec* Clone(int newnumvecs) const {
        return new MyMultiVec(length, newnumvecs);
    }

    MyMultiVec* CloneCopy() const {
        /* The newly created MyMultiVec will have one data block of size mylength*numvecs */
        MyMultiVec* v = new MyMultiVec(length, numvecs);

        /* Copy data blocks from this multi-vector into the single data block of v */
        int numblocks = (int) blocks.size();    // number of blocks in this multi-vector
        int vcount = 0;                         // number of vectors copied to v so far
        for(int b = 0; b < numblocks; b++) {
            int w = block_widths[b];
            int n = w * mylength;       // total size of block b
            const ScalarType* x = block_starts[b];
            ScalarType* y = v->block_starts[0] + mylength*vcount;
            blas.COPY(n, x, 1, y, 1);
            vcount += w;
        }

        return v;
    }

    MyMultiVec* CloneCopy(const std::vector<int>& index) const {
        /* The newly created MyMultiVec will have one data block of size mylength*newnumvecs */
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(length, newnumvecs);

        /* Copy individual vectors from this multi-vector into the single data block of v */
        for(int j = 0; j < newnumvecs; j++) {
            assert(index[j] >= 0 && index[j] < numvecs);

            /* Get info about vector number index[j] */
            VectorInfo info = index_to_block_map[index[j]];
            int b = info.b;     // block owning vector number index[j]
            int k = info.k;     // local index of this vector within block b

            /* Copy vector k of block b into v */
            int n = mylength;
            const ScalarType* x = block_starts[b] + mylength*k;
            ScalarType* y = v->block_starts[0] + mylength*j;
            blas.COPY(n, x, 1, y, 1);
        }

        return v;
    }

    const MyMultiVec* CloneView(const std::vector<int>& index) const {
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(length, 0);

        /* Build up v from the individual vectors of this multi-vector */
        int lastb = -1;         // last block (the one containing vector number index[j-1])
        int lastk = 0;          // last local index
        for(int j = 0; j < newnumvecs; j++) {
            assert(index[j] >= 0 && index[j] < numvecs);

            /* Get info about vector number index[j] */
            VectorInfo info = index_to_block_map[index[j]];
            int b = info.b;     // block owning vector number index[j]
            int k = info.k;     // local index of this vector within block b

            if(b == lastb && k == lastk + 1) {
                /* No need to add a new block, just widen the current one */
                v->block_widths.back() += 1;
            }
            else {
                /* Block not referenced before, or not contiguously; add new block */
                v->blocks.push_back(blocks[b]);
                v->block_widths.push_back(1);
                v->block_starts.push_back(block_starts[b] + mylength*k);
            }

            lastb = b;
            lastk = k;
        }

        /* Build vector-to-block mapping */
        for(int b = 0; b < (int) v->blocks.size(); b++) {
            int w = v->block_widths[b];
            for(int k = 0; k < w; k++) {
                VectorInfo info = { b, k };
                v->index_to_block_map.push_back(info);
            }
        }
        v->numvecs = newnumvecs;
        assert(v->numvecs == (int) v->index_to_block_map.size());
        return v;
    }

    MyMultiVec* CloneViewNonConst(const std::vector<int>& index) {
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(length, 0);

        /* Build up v from the individual vectors of this multi-vector */
        int lastb = -1;         // last block (the one containing vector number index[j-1])
        int lastk = 0;          // last local index
        for(int j = 0; j < newnumvecs; j++) {
            assert(index[j] >= 0 && index[j] < numvecs);

            /* Get info about vector number index[j] */
            VectorInfo info = index_to_block_map[index[j]];
            int b = info.b;     // block owning vector number index[j]
            int k = info.k;     // local index of this vector within block b

            if(b == lastb && k == lastk + 1) {
                /* No need to add a new block, just widen the current one */
                v->block_widths.back() += 1;
            }
            else {
                /* Block not referenced before, or not contiguously; add new block */
                v->blocks.push_back(blocks[b]);
                v->block_widths.push_back(1);
                v->block_starts.push_back(block_starts[b] + mylength*k);
            }

            lastb = b;
            lastk = k;
        }

        /* Build vector-to-block mapping */
        for(int b = 0; b < (int) v->blocks.size(); b++) {
            int w = v->block_widths[b];
            for(int k = 0; k < w; k++) {
                VectorInfo info = { b, k };
                v->index_to_block_map.push_back(info);
            }
        }
        v->numvecs = newnumvecs;
        assert(v->numvecs == (int) v->index_to_block_map.size());
        return v;
    }

    /* Attribute methods */

    int GetVecLength() const {
        /* Global vector length */
        return length;
    }

    int GetNumberVecs() const {
        return numvecs;
    }

    /* Update methods */

    void MvTimesMatAddMv(ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
            const Teuchos::SerialDenseMatrix<int,ScalarType>& B, ScalarType beta)
    {
        assert(length == A.GetVecLength());
        assert(B.numRows() == A.GetNumberVecs());
        assert(B.numCols() <= numvecs);

        int m = mylength;
        int k = A.GetNumberVecs();
        int n = numvecs;

        /* Cast A as a MyMultiVec (it is guaranteed to be one) */
        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL);
        assert(mylength == pA->mylength);

        /* In the most common usage case, both A and *this consist of a single
         * data block.  Optimize for this case. */
        if(pA->blocks.size() == 1 && this->blocks.size() == 1) {
            const ScalarType* a = pA->block_starts[0];
            const ScalarType* b = B.values();
            ScalarType* c = this->block_starts[0];

            /* BLAS calls can only be used if the memory regions associated
             * with a and c do not overlap. */
            if(!array_overlap(a, m*k, c, m*n)) {
                int lda = m;
                int ldb = B.stride();
                int ldc = m;
                blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
                return;
            }
        }

        /* Fall back to slow, row-by-row multiplication. */
        std::vector<ScalarType> u(k), v(n);
        for(int i = 0; i < m; i++) {
            /* Copy row i of multi-vector A into temporary u */
            for(int j = 0; j < k; j++)
                u[j] = pA->local(i,j);

            /* Copy row i of this multi-vector into temporary v */
            for(int j = 0; j < n; j++)
                v[j] = this->local(i,j);

            /* Compute vector-matrix product v = alpha*u*B + beta*v */
            const ScalarType* b = B.values();
            int ldb = B.stride();
            blas.GEMV(Teuchos::TRANS, k, n, alpha, b, ldb, &u[0], 1, beta, &v[0], 1);

            /* Copy temporary v back into this multi-vector */
            for(int j = 0; j < n; j++)
                this->local(i,j) = v[j];
        }
    }

    void MvAddMv(ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
                 ScalarType beta,  const Anasazi::MultiVec<ScalarType>& B)
    {
        assert(length == A.GetVecLength() && length == B.GetVecLength());
        assert(numvecs == A.GetNumberVecs() && numvecs == B.GetNumberVecs());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && mylength == pA->mylength);

        const MyMultiVec* pB = dynamic_cast<const MyMultiVec*>(&B);
        assert(pB != NULL && mylength == pB->mylength);

        /* Load B into C, scale by beta, then set C = alpha*A + C. */
        for(int j = 0; j < numvecs; j++) {
            const ScalarType* a = pA->vector(j);
            const ScalarType* b = pB->vector(j);
            ScalarType* c = this->vector(j);

//            for(int i = 0; i < mylength; i++)
//                c[i] = alpha*a[i] + beta*b[i];
            blas.COPY(mylength, b, 1, c, 1);
            blas.SCAL(mylength, beta, c, 1);
            blas.AXPY(mylength, alpha, a, 1, c, 1);
        }
    }

    void MvTransMv(ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A,
                   Teuchos::SerialDenseMatrix<int,ScalarType>& C) const
    {
        assert(length == A.GetVecLength());
        assert(numvecs <= C.numCols() && A.GetNumberVecs() <= C.numRows());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && mylength == pA->mylength);

        int m = pA->GetNumberVecs();
        int n = numvecs;
        int k = mylength;

        /* Local buffer for matrix elements */
        std::vector<ScalarType> buf(m*n);

        if(pA->blocks.size() == 1 && this->blocks.size() == 1) {
            /* Common case: each multi-vector consists of a single block. */
            const ScalarType* a = pA->block_starts[0];
            const ScalarType* b = this->block_starts[0];
            ScalarType* c = &buf[0];
            int lda = mylength;
            int ldb = mylength;
            int ldc = m;
            blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, k, 1, a, lda, b, ldb, 0, c, ldc);
        }
        else {
            /* General case: multiply column-by-column */
            for(int i = 0; i < m; i++) {
                const ScalarType* a = pA->vector(i);
                for(int j = 0; j < n; j++) {
                    const ScalarType* b = this->vector(j);
                    buf[i + m*j] = blas.DOT(k, a, 1, b, 1);
                }
            }
        }

#ifdef OPSEC_USE_MPI
        /* Sum contributions from all processes */
        MPI_Allreduce(MPI_IN_PLACE, &buf[0], m*n,
                      Teuchos::RawMPITraits<ScalarType>::type(), MPI_SUM, MPI_COMM_WORLD);
#endif

        /* Copy local buffer back into C */
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                C(i,j) = buf[i + m*j];
    }

    void MvDot(const Anasazi::MultiVec<ScalarType>& A, std::vector<ScalarType>& c) const {
        assert(length == A.GetVecLength() && numvecs == A.GetNumberVecs());
        assert(numvecs <= (int) c.size());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && mylength == pA->mylength);

        int n = mylength;
        for(int j = 0; j < numvecs; j++) {
            const ScalarType* a = pA->vector(j);
            const ScalarType* b = this->vector(j);
            c[j] = blas.DOT(n, a, 1, b, 1);
        }

#ifdef OPSEC_USE_MPI
        /* Sum contributions from all processes */
        MPI_Allreduce(MPI_IN_PLACE, &c[0], numvecs,
                      Teuchos::RawMPITraits<ScalarType>::type(), MPI_SUM, MPI_COMM_WORLD);
#endif
    }


    /* Norm */

    void MvNorm(std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec) const
    {
        typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
        assert(numvecs <= (int) normvec.size());

        int n = mylength;
        /* Compute local sum of squares */
        for(int j = 0; j < numvecs; j++) {
            const ScalarType* x = this->vector(j);
            normvec[j] = blas.DOT(n, x, 1, x, 1);
        }

#ifdef OPSEC_USE_MPI
        /* Sum contributions from all processes */
        MPI_Allreduce(MPI_IN_PLACE, &normvec[0], numvecs,
                      Teuchos::RawMPITraits<MagnitudeType>::type(), MPI_SUM, MPI_COMM_WORLD);
#endif

        /* Take square root */
        for(int j = 0; j < numvecs; j++)
            normvec[j] = Teuchos::ScalarTraits<MagnitudeType>::squareroot(normvec[j]);
    }


    /* Initialization methods */

    void SetBlock(const Anasazi::MultiVec<ScalarType>& A, const std::vector<int>& index) {
        assert(length == A.GetVecLength() && A.GetNumberVecs() >= (int) index.size());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && mylength == pA->mylength);

        int n = mylength;
        for(int j = 0; j < (int) index.size(); j++) {
            const ScalarType* a = pA->vector(j);
            ScalarType* b = this->vector(index[j]);
            blas.COPY(n, a, 1, b, 1);
        }
    }

    void MvScale(ScalarType alpha) {
        /* Scale elements block-by-block */
        for(int b = 0; b < (int) blocks.size(); b++) {
            int w = block_widths[b];
            int n = w*mylength;
            ScalarType* x = block_starts[b];
            blas.SCAL(n, alpha, x, 1);
        }
    }

    void MvScale(const std::vector<ScalarType>& alpha) {
        assert(numvecs <= (int) alpha.size());
        int n = mylength;
        for(int j = 0; j < numvecs; j++) {
            ScalarType* x = this->vector(j);
            blas.SCAL(n, alpha[j], x, 1);
        }
    }

    void MvRandom() {
        for(int j = 0; j < numvecs; j++) {
            ScalarType* x = this->vector(j);
            for(int i = 0; i < mylength; i++)
                x[i] = Teuchos::ScalarTraits<ScalarType>::random();
        }
    }

    void MvInit(ScalarType alpha) {
        for(int j = 0; j < numvecs; j++) {
            ScalarType* x = this->vector(j);
            for(int i = 0; i < mylength; i++)
                x[i] = alpha;
        }
    }

    void MvPrint(std::ostream& os) const {
#ifdef OPSEC_USE_MPI
        MPI_Datatype mpitype = Teuchos::RawMPITraits<ScalarType>::type();
        MPI_Status status;
        int goahead = 1;
#endif

        /* Copy local data into contiguous buffer */
        ScalarType* buf = (ScalarType*) malloc(mylength*numvecs*sizeof(ScalarType));
        for(int i = 0; i < mylength; i++)
            for(int j = 0; j < numvecs; j++)
                buf[j*mylength + i] = local(i,j);

        if(me == 0) {
            /* Root process: gather data from each process and print */

            os << "Object MyMultiVec" << std::endl;
            os << "  length = " << length << ", numvecs = " << numvecs << std::endl;

            for(int p = 0; p < nprocs; p++) {
                /* Local length for process p (always <= mylength for root process) */
                int plength = (length/nprocs) + (p < (length % nprocs));

#ifdef OPSEC_USE_MPI
                if(p != 0) {
                    /* Signal process p that we're ready for its data, and wait for response */
                    MPI_Send(&goahead, 1, MPI_INT, p, 0, MPI_COMM_WORLD);
                    MPI_Recv(buf, plength*numvecs, mpitype, p, 0, MPI_COMM_WORLD, &status);
                }
#endif

                /* Print data */
                os << "  Process " << p << ":" << std::endl;
                for(int i = 0; i < plength; i++) {
                    os << "    ";
                    for(int j = 0; j < numvecs; j++)
                        os << buf[j*plength + i] << " ";
                    os << std::endl;
                }
            }
        }
#ifdef OPSEC_USE_MPI
        else {
            /* Non-root process: wait for signal, then send local data */
            MPI_Recv(&goahead, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send(buf, mylength*numvecs, mpitype, 0, 0, MPI_COMM_WORLD);
        }
#endif

        free(buf);
    }


    /***** Additional methods *****/

    int GetLocalLength() const {
        /* Local vector length */
        return mylength;
    }

    int GetLocalStart() const {
        return mystart;
    }


/* For testing and debugging */
//private:
public:

    int length;
    int numvecs;
    int nprocs;         // number of processes in MPI communicator
    int me;             // my process rank
    int mylength;       // number of elements assigned to me
    int mystart;        // index of first element assigned to me

    Teuchos::BLAS<int, ScalarType> blas;        // BLAS interface

    /* A data block is a contiguous region of memory, representing m vectors
     * each of length n.  The data is arranged in memory such that each
     * length-n vector is contiguous (i.e. it is a column-major n-by-m matrix).
     * A MultiVec is made up of one or more data blocks.  Each DataBlock is
     * reference-counted to facilitate creating views from one MultiVec into
     * another. */
    class DataBlock {
        /* Internal data structure used by DataBlock */
        struct RawDataBlock {
            ScalarType* p;
            int refcount;
        };

        int n;                  // length of each vector in memory
        int m;                  // number of vectors in block
        RawDataBlock* raw;      // reference-counted pointer to raw data

    public:
        DataBlock() {
            n = m = 0;
            raw = new RawDataBlock;
            raw->p = NULL;
            raw->refcount = 1;
        }

        DataBlock(int n_, int m_) {
            n = n_;
            m = m_;
            assert(n >= 0 && m >= 0);
            raw = new RawDataBlock;
            raw->p = (ScalarType*) malloc(n*m*sizeof(ScalarType));
            raw->refcount = 1;
        }

        DataBlock(const DataBlock& b) {
            n = b.n;
            m = b.m;
            raw = b.raw;
            ++raw->refcount;
        }

        ~DataBlock() {
            if(--raw->refcount <= 0) {
                free(raw->p);
                delete raw;
            }
        }

        DataBlock& operator=(const DataBlock& b) {
            n = b.n;
            m = b.m;
            if(raw != b.raw) {  // avoid problems with self-assignment
                if(--raw->refcount <= 0) {
                    free(raw->p);
                    delete raw;
                }
                raw = b.raw;
                ++raw->refcount;
            }
            return *this;
        }

        operator ScalarType*()             { return raw->p; }
        operator const ScalarType*() const { return raw->p; }
    };

    /* Array of data blocks making up this MultiVec */
    std::vector<DataBlock> blocks;

    /* Number of vectors used in each data block (note that block_widths[i] may
     * differ from blocks[i].m, e.g. if the current MultiVec is only a partial
     * view into another MultiVec, because the entire data block might
     * not be used by the current MultiVec). */
    std::vector<int> block_widths;

    /* Pointer to first used element of each data block */
    std::vector<ScalarType*> block_starts;

    struct VectorInfo {
        int b;  // block number
        int k;  // local index of vector within block
    };

    /* Mapping that gives the block number b, and intra-block index k,
     * corresponding to a given vector */
    std::vector<VectorInfo> index_to_block_map;
};

#endif // MYMULTIVEC_H
