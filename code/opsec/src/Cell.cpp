#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <stdio.h>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "Cell.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"

Cell* ReadCells(const char* cellfile, int *Ncells_) {
    int Ncells;
    Cell* cells;

    int nprocs = 1, me = 0;
#ifdef OPSEC_USE_MPI
    int ready;
    MPI_Initialized(&ready);
    if(ready) {
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
    }
#endif

    /* Have root process read in cells from file */
    if(me == 0) {
        size_t n, size;
        char endian, fmt[ABN_MAX_FORMAT_LENGTH];
        Config cellopts = cfg_new();

        /* Read basis cells from file */
        FILE* fp = fopen(cellfile, "r");
        if(fp == NULL) {
            opsec_error("ReadCells: could not open '%s'\n", cellfile);
            opsec_abort(1);
        }
        if(abn_read(fp, (void**) &cells, &n, &size, &endian, fmt, cellopts) != 0) {
            opsec_error("ReadCells: error reading cells from '%s'\n", cellfile);
            opsec_abort(1);
        }

        Ncells = cfg_get_int(cellopts, "Ncells");
        if((int)n != Ncells) {
            opsec_error("ReadCells: consistency check fail for Ncells: %d != %zd\n", Ncells, n);
            opsec_abort(1);
        }
        if(size != sizeof(Cell)) {
            opsec_error("ReadCells: consistency check fail for sizeof(Cell): %zd != %zd\n", sizeof(Cell), size);
            opsec_abort(1);
        }

        opsec_debug("Read %d basis cells from '%s'\n", Ncells, cellfile);

        fclose(fp);
        cfg_destroy(cellopts);
    }

#ifdef OPSEC_USE_MPI
    if(nprocs > 1) {
        /* Broadcast cell data to all other processes */
        MPI_Bcast(&Ncells, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(me != 0) {
            cells = (Cell*) opsec_malloc(Ncells*sizeof(Cell));
            if(cells == NULL) {
                opsec_error("ReadCells: could not allocate memory for cells on process %d\n", me);
                opsec_abort(1);
            }
        }
        MPI_Bcast(cells, Ncells*sizeof(Cell), MPI_BYTE, 0, MPI_COMM_WORLD);
    }
#endif

    if(Ncells_) *Ncells_ = Ncells;
    return cells;
}
