#include <cassert>
#include <vector>

#include <AnasaziMultiVec.hpp>
#include <AnasaziOperator.hpp>

#include "opsec.h"
#include "slp.h"

/* Return n divided by p, rounded up to nearest integer */
static inline int divup(int n, int p) {
    return (n + p - 1)/p;
}

/****** ColumnBlock ***********************************************************
 * A utility class representing a contiguous span of columns of a distributed
 * matrix.  The matrix distribution is presumed to be one-dimensional, i.e.
 * each process is given full rows of the matrix.
 *
 * The example below shows a matrix with 7 columns, broken into three blocks.
 * The first block spans columns 0 and 1, the second blocks spans columns 2-5,
 * and the third block contains column 7.  The parameters width and start for
 * each block indicate the positioning of the block within the matrix.
 *    0 1 2 3 4 5 6     0 1 2 3 4 5 6 
 *   ---------------   ---------------
 *   |             |   |   :       : |
 *   |     raw     |   |   :   X   : |
 *   |    block    |   |   :       : |
 *   |             |   |   :       : |
 *   ---------------   ---------------
 *                     block X within raw block
 *                      width = 4, start = 2
 */
struct ColumnBlock {
    /* Internal data structure used by ColumnBlock.  A raw block is an actual
     * block of local memory (values) and the array descriptor that goes
     * with it.  Raw blocks are reference-counted to allow sharing. */
    struct RawBlock {
        int refcount;
        slp::Descriptor* desc;
        real* values;
    };

    int width;      // number of columns in this block
    int start;      // starting column of this block within the raw block
    int lld;        // memory stride between columns in local memory
    RawBlock* raw;  // raw block this block points into


    /* Default constructor: Create an empty block. */
    ColumnBlock() {
        width = 0;
        start = 0;
        raw = new RawBlock;
        raw->refcount = 1;
        raw->desc = NULL;
        raw->values = NULL;
        lld = 0;
    }

    /* Create a block for an n-by-k matrix, with row-wise blocking factor nb. */
    ColumnBlock(const slp::Context* pcontext, int n, int k, int nb) {
        assert(n > 0 && k > 0 && nb > 0);
        width = k;
        start = 0;
        raw = new RawBlock;
        raw->refcount = 1;
        raw->desc = pcontext->new_descriptor(n, k, nb, k);
        raw->values = (real*) malloc(raw->desc->local_size() * sizeof(real));
        lld = raw->desc->lld;
    }

    ColumnBlock(const ColumnBlock& b) {
        width = b.width;
        start = b.start;
        lld = b.lld;
        raw = b.raw;
        ++raw->refcount;
    }

    ~ColumnBlock() {
        if(--raw->refcount <= 0) {
            free(raw->values);
            delete raw->desc;
            delete raw;
        }
    }

    ColumnBlock& operator=(const ColumnBlock& b) {
        width = b.width;
        start = b.start;
        lld = b.lld;
        if(raw != b.raw) {  // avoid problems with self-assignment
            if(--raw->refcount <= 0) {
                free(raw->values);
                delete raw->desc;
                delete raw;
            }
            raw = b.raw;
            ++raw->refcount;
        }
        return *this;
    }

//    operator real*()             { return raw->values; }
//    operator const real*() const { return raw->values; }

    /* Pointer to column d within this block (0 <= d < width).
     * Returns NULL if block is uninitialized. */
    real* column(int d) {
        assert(0 <= d && d < width);
        return raw->values + lld*(start + d);
    }
    const real* column(int d) const {
        assert(0 <= d && d < width);
        return raw->values + lld*(start + d);
    }

    /* Return the slp::Matrix that this block belongs to. */
    slp::Matrix<real> matrix() const { return slp::Matrix<real>(raw->desc, raw->values); }
};


/***** MyMultiVec *************************************************************
 * A distributed-memory implementation of Anasazi's MultiVec interface, using
 * ScaLAPACK. */
class MyMultiVec : public Anasazi::MultiVec<real> {
public:
    /* Constructor
     *   pcontext = BLACS process grid
     *   length = length of vectors
     *   numvecs = number of vectors
     *   nb = blocking factor for length
     */
    MyMultiVec(const slp::Context* pcontext, int length, int numvecs = 0, int nb = 0)
        : pcontext(pcontext), length(length), numvecs(numvecs), nb(nb)
    {
        assert(pcontext != NULL);
        assert(length > 0);

        /* Choose row blocking factor */
        if(nb <= 0)
            nb = divup(length, pcontext->nprow);

        /* Use one block for the entire multi-vector */
        if(numvecs > 0) {
            blocks.push_back(ColumnBlock(pcontext, length, numvecs, nb));
            for(int j = 0; j < numvecs; j++)
                column_map.push_back(ColumnPointer(0, j));
        }
    }

    /* Destructor */
    ~MyMultiVec() {
        /* Cleanup should be automatic */
    }


    /***** Accessors **********************************************************/

    /* Access local vectors */
    real* vector(int j) {
        assert(0 <= j && j < GetNumberVecs());
        ColumnPointer cp = column_map[j];
        return blocks[cp.b].column(cp.d);
    }
    const real* vector(int j) const {
        assert(0 <= j && j < GetNumberVecs());
        ColumnPointer cp = column_map[j];
        return blocks[cp.b].column(cp.d);
    }

    /* Access elements by local row index iloc */
    real& local(int iloc, int j) {
        assert(0 <= iloc && iloc < GetLocalLength());
        return vector(j)[iloc];
    }
    const real& local(int iloc, int j) const {
        assert(0 <= iloc && iloc < GetLocalLength());
        return vector(j)[iloc];
    }


    /***** Anasazi::MultiVec<real> interface **********************************/

    /* Creation methods */

    MyMultiVec* Clone(int newnumvecs) const {
        return new MyMultiVec(pcontext, length, newnumvecs, nb);
    }

    MyMultiVec* CloneCopy() const {
        /* The newly created MyMultiVec will have one data block of width numvecs */
        MyMultiVec* v = new MyMultiVec(pcontext, length, numvecs, nb);
        slp::Matrix<real> vmat = v->blocks[0].matrix(); // matrix defining single block of v
        int vcount = 0; // number of vectors copied to v so far

        /* Copy data blocks from this multi-vector into the single block of v */
        for(std::vector<ColumnBlock>::const_iterator bi = blocks.begin(); bi != blocks.end(); ++bi) {
            slp::redistribute(length, bi->width,
                              bi->matrix(), 0, bi->start,
                              vmat, 0, vcount);
            vcount += bi->width;
        }

        assert(vcount == numvecs);
        return v;
    }

    MyMultiVec* CloneCopy(const std::vector<int>& index) const {
        /* The newly created MyMultiVec will have one data block of width newnumvecs */
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(pcontext, length, newnumvecs, nb);
        slp::Matrix<real> vmat = v->blocks[0].matrix();

        /* Copy individual vectors from this multi-vector into the single data block of v */
        for(int j = 0; j < newnumvecs; j++) {
            assert(0 <= index[j] && index[j] < numvecs);

            /* Get info about column number index[j] */
            ColumnPointer cp = column_map[index[j]];
            int b = cp.b;       // index of block owning column index[j]
            int d = cp.d;       // local index of this column within block b

            /* Copy column d of block b into v (as column j of block 0) */
            slp::redistribute(length, 1,
                              blocks[b].matrix(), 0, d + blocks[b].start,
                              vmat, 0, j);
        }

        return v;
    }

    const MyMultiVec* CloneView(const std::vector<int>& index) const {
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(pcontext, length, 0, nb);
                        // start with 0 vectors, build up one at a time

        /* Build up v from the individual vectors of this multi-vector */
        int lastb = -1;         // last block (the one containing vector number index[j-1])
        int lastd = 0;          // last local index
        for(int j = 0; j < newnumvecs; j++) {
            assert(0 <= index[j] && index[j] < numvecs);

            /* Find column number index[j] */
            ColumnPointer cp = column_map[index[j]];
            int b = cp.b;       // block owning column index[j]
            int d = cp.d;       // local index of this column within block b

            if(b == lastb && d == lastd + 1) {
                /* No need to add a new block, just widen the current one */
                v->blocks.back().width += 1;
            }
            else {
                /* ColumnBlock not referenced before, or not contiguously; add new block */
                ColumnBlock newblock = blocks[b];
                newblock.width = 1;
                newblock.start = d;
                v->blocks.push_back(newblock);
            }

            lastb = b;
            lastd = d;
        }

        /* Build column-to-block mapping from scratch */
        int numnewblocks = (int) v->blocks.size();
        for(int b = 0; b < numnewblocks; b++) {
            int w = v->blocks[b].width;
            for(int d = 0; d < w; d++)
                v->column_map.push_back(ColumnPointer(b, d));
        }
        v->numvecs = newnumvecs;
        assert(v->numvecs == (int) v->column_map.size());
        return v;
    }

    MyMultiVec* CloneViewNonConst(const std::vector<int>& index) {
        int newnumvecs = (int) index.size();
        MyMultiVec* v = new MyMultiVec(pcontext, length, 0, nb);
                        // start with 0 vectors, build up one at a time

        /* Build up v from the individual vectors of this multi-vector */
        int lastb = -1;         // last block (the one containing vector number index[j-1])
        int lastd = 0;          // last local index
        for(int j = 0; j < newnumvecs; j++) {
            assert(0 <= index[j] && index[j] < numvecs);

            /* Find column number index[j] */
            ColumnPointer cp = column_map[index[j]];
            int b = cp.b;       // block owning column index[j]
            int d = cp.d;       // local index of this column within block b

            if(b == lastb && d == lastd + 1) {
                /* No need to add a new block, just widen the current one */
                v->blocks.back().width += 1;
            }
            else {
                /* ColumnBlock not referenced before, or not contiguously; add new block */
                ColumnBlock newblock = blocks[b];
                newblock.width = 1;
                newblock.start = d;
                v->blocks.push_back(newblock);
            }

            lastb = b;
            lastd = d;
        }

        /* Build column-to-block mapping from scratch */
        int numnewblocks = (int) v->blocks.size();
        for(int b = 0; b < numnewblocks; b++) {
            int w = v->blocks[b].width;
            for(int d = 0; d < w; d++)
                v->column_map.push_back(ColumnPointer(b, d));
        }
        v->numvecs = newnumvecs;
        assert(v->numvecs == (int) v->column_map.size());
        return v;
    }

    /* Attribute methods */

    int GetVecLength() const {
        return length;
    }

    int GetNumberVecs() const {
        return numvecs;
    }

    int GetLocalLength() const {
        return pcontext->numrows(length, nb);
    }


    /* Update methods */

    /* C = alpha*A*B + beta*C   (where C = *this) */
    void MvTimesMatAddMv(real alpha, const Anasazi::MultiVec<real>& A,
                         const Teuchos::SerialDenseMatrix<int,real>& B, real beta)
    {
        assert(GetVecLength() == A.GetVecLength());
        assert(B.numRows() == A.GetNumberVecs());
        assert(B.numCols() <= GetNumberVecs());

        int m = GetLocalLength();
        int n = A.GetNumberVecs();
        int k = B.numCols();

        /* Cast A as a MyMultiVec (it is guaranteed to be one) */
        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && pA->GetLocalLength() == GetLocalLength());

        /* In the most common usage case, both A and C consist of a single
         * block.  Optimize for this case. */
        if(pA->blocks.size() == 1 && this->blocks.size() == 1) {
            const ColumnBlock& Ablock = pA->blocks.front();
            ColumnBlock& Cblock = this->blocks.front();

            const real* a = Ablock.column(0);
            const real* b = B.values();
            real* c = Cblock.column(0);

            int lda = Ablock.lld;
            int ldb = B.stride();
            int ldc = Cblock.lld;

            blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        }
        else {
            /* Fall back to (slow!) element-wise multiplication. */
            for(int i = 0; i < m; i++) {
                for(int j = 0; j < k; j++) {
                    double tmp = 0;
                    for(int p = 0; p < n; p++)
                        tmp += pA->local(i,p) * B(p,j);
                    this->local(i,j) = alpha*tmp + beta*this->local(i,j);
                }
            }
        }
    }

    /* C = alpha*A + beta*B  (where C = *this) */
    void MvAddMv(real alpha, const Anasazi::MultiVec<real>& A,
                 real beta,  const Anasazi::MultiVec<real>& B)
    {
        assert(GetVecLength() == A.GetVecLength() && GetVecLength() == B.GetVecLength());
        assert(GetNumberVecs() == A.GetNumberVecs() && GetNumberVecs() == B.GetNumberVecs());

        /* Cast to my types */
        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        const MyMultiVec* pB = dynamic_cast<const MyMultiVec*>(&B);
        assert(pA != NULL && pB != NULL);
        assert(GetLocalLength() == pA->GetLocalLength() && GetLocalLength() == pB->GetLocalLength());

        int m = GetLocalLength();
        int n = GetNumberVecs();

        /* Perform the addition one column at a time */
        for(int j = 0; j < n; j++) {
            const real* a = pA->vector(j);
            const real* b = pB->vector(j);
            real* c = this->vector(j);

            for(int i = 0; i < m; i++)
                c[i] = alpha*a[i] + beta*b[i];
//            blas.COPY(nloc, b, 1, c, 1);
//            blas.SCAL(nloc, beta, c, 1);
//            blas.AXPY(nloc, alpha, a, 1, c, 1);
        }
    }

    /* C = alpha*A^T*B   (where B = *this) */
    void MvTransMv(real alpha, const Anasazi::MultiVec<real>& A,
                   Teuchos::SerialDenseMatrix<int,real>& C) const
    {
        assert(GetVecLength() == A.GetVecLength());
        assert(GetNumberVecs() <= C.numCols() && A.GetNumberVecs() <= C.numRows());

        /* Cast to my type */
        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL);
        assert(GetLocalLength() == pA->GetLocalLength());
        assert(pcontext == pA->pcontext);

        int m = A.GetNumberVecs();
        int n = this->GetLocalLength();
        int k = this->GetNumberVecs();
        int ldc = C.stride();

        if(pA->blocks.size() == 1 && this->blocks.size() == 1) {
            /* Common case: each multi-vector consists of a single block,
             * perform the local multiplication in one step */
            const ColumnBlock& Ablock = pA->blocks[0];
            const ColumnBlock& Bblock = this->blocks[0];

            const real* a = Ablock.column(0);
            const real* b = Bblock.column(0);
            real* c = C.values();
            int lda = Ablock.lld;
            int ldb = Bblock.lld;
            blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, m, n, k, alpha, a, lda, b, ldb, 0.0, c, ldc);
        }
        else {
            /* General case: multiply column-by-column */
            for(int i = 0; i < m; i++) {
                const real* a = pA->vector(i);
                for(int j = 0; j < k; j++) {
                    const real* b = this->vector(j);
                    C(i,j) = blas.DOT(n, a, 1, b, 1);
                }
            }
        }

        /* Sum contributions from all processes (all-reduce) */
        slp::gsum2d(pcontext, "All", " ", m, k, C.values(), ldc, -1, -1);
    }

    void MvDot(const Anasazi::MultiVec<real>& A, std::vector<real>& c) const {
        assert(GetVecLength() == A.GetVecLength() && GetNumberVecs() == A.GetNumberVecs());
        assert(GetNumberVecs() <= (int) c.size());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && GetLocalLength() == pA->GetLocalLength() && pcontext == pA->pcontext);

        int m = GetLocalLength();
        int n = GetNumberVecs();
        for(int j = 0; j < n; j++) {
            const real* a = pA->vector(j);
            const real* b = this->vector(j);
            c[j] = blas.DOT(m, a, 1, b, 1);
        }

        /* Sum contributions from all processes (broadcast to all processes) */
        slp::gsum2d(pcontext, "All", " ", n, 1, &c[0], 1, -1, -1);
    }


    /* Norm */

    void MvNorm(std::vector<real>& normvec) const {
        assert(GetNumberVecs() <= (int) normvec.size());

        int m = GetLocalLength();
        int n = GetNumberVecs();

        /* Compute local sum of squares */
        for(int j = 0; j < n; j++) {
            const real* x = this->vector(j);
            normvec[j] = blas.DOT(m, x, 1, x, 1);
        }

        /* Sum contributions from all processes (broadcast to all processes) */
        slp::gsum2d(pcontext, "All", " ", n, 1, &normvec[0], 1, -1, -1);

        /* Take square root */
        for(int j = 0; j < n; j++)
            normvec[j] = sqrt(normvec[j]);
    }


    /* Initialization methods */

    void SetBlock(const Anasazi::MultiVec<real>& A, const std::vector<int>& index) {
        assert(GetVecLength() == A.GetVecLength() && A.GetNumberVecs() >= (int) index.size());

        const MyMultiVec* pA = dynamic_cast<const MyMultiVec*>(&A);
        assert(pA != NULL && GetLocalLength() == pA->GetLocalLength());

        int m = GetLocalLength();
        int n = (int) index.size();
        for(int j = 0; j < n; j++) {
            const real* a = pA->vector(j);
            real* b = this->vector(index[j]);
            blas.COPY(m, a, 1, b, 1);
        }
    }

    void MvScale(real alpha) {
        int m = GetLocalLength();
        int n = GetNumberVecs();
        for(int j = 0; j < n; j++) {
            real* x = this->vector(j);
            blas.SCAL(m, alpha, x, 1);
        }
    }

    void MvScale(const std::vector<real>& alpha) {
        assert(GetNumberVecs() <= (int) alpha.size());
        int m = GetLocalLength();
        int n = GetNumberVecs();
        for(int j = 0; j < n; j++) {
            real* x = this->vector(j);
            blas.SCAL(m, alpha[j], x, 1);
        }
    }

    void MvRandom() {
        int m = GetLocalLength();
        int n = GetNumberVecs();
        for(int j = 0; j < n; j++) {
            real* x = this->vector(j);
            for(int i = 0; i < m; i++)
                x[i] = Teuchos::ScalarTraits<real>::random();
        }
    }

    void MvInit(real alpha) {
        int m = GetLocalLength();
        int n = GetNumberVecs();
        for(int j = 0; j < n; j++) {
            real* x = this->vector(j);
            for(int i = 0; i < m; i++)
                x[i] = alpha;
        }
    }

    void MvPrint(std::ostream& os) const {
        os << "MyMultiVec::MvPrint() not implemented yet." << std::endl;
#if 0
#ifdef OPSEC_USE_MPI
        MPI_Datatype mpitype = Teuchos::RawMPITraits<real>::type();
        MPI_Status status;
        int goahead = 1;
#endif

        /* Copy local data into contiguous buffer */
        real* buf = (real*) malloc(mylength*numvecs*sizeof(real));
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
#endif
    }


/* For testing and debugging */
//private:
public:
    const slp::Context* pcontext;       // BLACS process grid
    int length;                         // number of rows in each vector
    int numvecs;                        // number of multi-vectors
    int nb;                             // row-wise blocking factor

    Teuchos::BLAS<int, real> blas;      // BLAS interface

    /* List of data blocks making up this MyMultiVec */
    std::vector<ColumnBlock> blocks;

    struct ColumnPointer {
        int b;  // index of block owning column
        int d;  // local index of column within block b
        ColumnPointer(int b = 0, int d = 0) : b(b), d(d) {}
    };

    /* Mapping that gives the block number b, and intra-block index j,
     * corresponding to a given column */
    std::vector<ColumnPointer> column_map;

    friend class MyOperator;
};


class MyOperator : public Anasazi::Operator<real> {
public:
    /* n is the global size of the matrix, pcontext is the BLACS process grid,
     * mb and nb are the row-wise and column-wise blocking factors */
    MyOperator(const slp::Context* pcontext, int n, MatrixFactory& matfact, int mb = 0, int nb = 0)
        : pcontext(pcontext), n(n)
    {
        /* Choose default blocking factors */
        if(mb <= 0)
            mb = divup(n, pcontext->nprow);
        if(nb <= 0)
            nb = divup(n, pcontext->npcol);

        /* Allocate memory for distributed n-by-n matrix */
        Adesc = pcontext->new_descriptor(n, n, mb, nb);
        Avalues = (real*) opsec_malloc(Adesc->local_size() * sizeof(real));

        /* Fill in local matrix values */
        int nrows = Adesc->num_local_rows();
        int ncols = Adesc->num_local_cols();
        std::vector<int> rows(nrows), cols(ncols);
        for(int i = 0; i < nrows; i++)
            rows[i] = Adesc->row_l2g(i);
        for(int j = 0; j < ncols; j++)
            cols[j] = Adesc->col_l2g(j);
        matfact.ComputeMatrixValues(rows, cols, Avalues, Adesc->lld);
    }

    /* Destructor */
    ~MyOperator() {
        free(Avalues);
        delete Adesc;
    }

    void Apply(const Anasazi::MultiVec<real>& x, Anasazi::MultiVec<real>& y) const {
        assert(x.GetVecLength() == n && y.GetVecLength() == n);
        assert(x.GetNumberVecs() == y.GetNumberVecs());

        /* Cast to my types */
        const MyMultiVec* px = dynamic_cast<const MyMultiVec*>(&x);
        MyMultiVec* py = dynamic_cast<MyMultiVec*>(&y);
        assert(px != NULL && py != NULL);

        slp::Matrix<real> A(Adesc, Avalues);

        int numvecs = px->GetNumberVecs();
        if(numvecs <= 0)
            return;

        /* Iterate over the blocks of each multi-vector.  Identify common
         * sub-blocks, and perform the multiplication in chunks.
         *   x   | . . . | . | | . . . |   blocks of x
         *
         *   y   | | . . . | | . . | . |   blocks of y
         *
         *   ->  | | . . | | | | . | . |   common sub-blocks of x and y
         */
        std::vector<ColumnBlock>::const_iterator bx = px->blocks.begin();
        std::vector<ColumnBlock>::iterator by = py->blocks.begin();
        int dx = 0;     // current starting position within block *bx
        int dy = 0;     // current starting position within block *by
        int w = 1;      // width of current common block
        for(int j = 0; j < numvecs; j++) {
            if(dx + w == bx->width || dy + w == by->width) {
                /* We're at the boundary of a common block.  So multiply A by
                 * the n-by-w sub-block of x, storing the result in the
                 * corresponding sub-block of y. */
                slp::Matrix<real> xmat = bx->matrix();
                slp::Matrix<real> ymat = by->matrix();
                slp::multiply(A, xmat, ymat,
                              'N', 'N',
                              n, w, n,
                              1.0, 0.0,
                              0, 0, 
                              0, dx + bx->start,
                              0, dy + by->start);
                /* Move the indices forward */
                dx += w;
                dy += w;
                /* Reset whichever boundary we hit */
                if(dx == bx->width) {
                    ++bx;
                    dx = 0;
                }
                if(dy == by->width) {
                    ++by;
                    dy = 0;
                }
                w = 1;
            }
            else {
                /* We're not at a boundary, so extend the common sub-block. */
                ++w;
            }
        }
        assert(bx == px->blocks.end() && by == py->blocks.end());
    }

private:
    /* BLACS process grid */
    const slp::Context* pcontext;

    /* Global matrix size */
    int n;

    /* Descriptor for distributed n-by-n matrix */
    slp::Descriptor* Adesc;

    /* Local values of distributed n-by-n matrix */
    real* Avalues;
};
