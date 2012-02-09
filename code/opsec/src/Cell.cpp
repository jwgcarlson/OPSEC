#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif

#include "Cell.h"
#include "abn.h"
#include "cfg.h"
#include "opsec.h"

#if 0
Cell* ReadCells(const char* cellfile, int* Ncells_) {
    size_t n, size;
    char endian, fmt[ABN_MAX_FORMAT_LENGTH];
    int Ncells;
    FILE* fp;
    Cell* cells = NULL;
    Config cellopts = cfg_new();

    fp = fopen(cellfile, "r");
    if(fp == NULL) {
        fprintf(stderr, "ReadCells: could not read cells from '%s'\n", cellfile);
        fflush(stderr);
        return NULL;
    }
    if(abn_read(fp, (void**) &cells, &n, &size, &endian, fmt, cellopts) != 0) {
        fprintf(stderr, "ReadCells: error reading cells from '%s'\n", cellfile);
        fflush(stderr);
        return NULL;
    }

    Ncells = cfg_get_int(cellopts, "Ncells");
    if((int)n != Ncells) {
        fprintf(stderr, "ReadCells: consistency check fail for Ncells: %d != %zd\n", Ncells, n);
        fflush(stderr);
    }
    if(size != sizeof(Cell)) {
        fprintf(stderr, "ReadCells: consistency check fail for sizeof(Cell): %zd != %zd\n", sizeof(Cell), size);
        fflush(stderr);
    }

    fclose(fp);
    cfg_destroy(cellopts);
    if(Ncells_) *Ncells_ = Ncells;
    return cells;
}
#endif

Cell* ReadCells(const char* cellfile, int *Ncells_) {
    int Ncells;
    Cell* cells;

    int nprocs = 1, me = 0, usempi = 0;
#ifdef OPSEC_USE_MPI
    MPI_Initialized(&usempi);
    if(usempi) {
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
            fprintf(stderr, "ReadCells: could not open '%s'\n", cellfile);
            fflush(stderr);
            opsec_abort(1);
        }
        if(abn_read(fp, (void**) &cells, &n, &size, &endian, fmt, cellopts) != 0) {
            fprintf(stderr, "ReadCells: error reading cells from '%s'\n", cellfile);
            fflush(stderr);
            opsec_abort(1);
        }

        Ncells = cfg_get_int(cellopts, "Ncells");
        if((int)n != Ncells) {
            fprintf(stderr, "ReadCells: consistency check fail for Ncells: %d != %zd\n", Ncells, n);
            fflush(stderr);
        }
        if(size != sizeof(Cell)) {
            fprintf(stderr, "ReadCells: consistency check fail for sizeof(Cell): %zd != %zd\n", sizeof(Cell), size);
            fflush(stderr);
        }

        printf("(TRACE) Read %d basis cells from '%s'\n", Ncells, cellfile); fflush(stdout);

        fclose(fp);
        cfg_destroy(cellopts);
    }


#ifdef OPSEC_USE_MPI
    if(usempi) {
        printf("(TRACE) Sending cells to other processes...\n"); fflush(stdout);

        /* Define MPI datatype for cells (this must be updated if struct Cell changes!) */
        MPI_Datatype cell_datatype;
        int blocklengths[2] = { 2, 8 };
        MPI_Aint displacements[2] = { 0, 8 };
        MPI_Datatype types[2] = { MPI_INT, MPI_DOUBLE };
        MPI_Type_create_struct(2, blocklengths, displacements, types, &cell_datatype);
        MPI_Type_commit(&cell_datatype);

        /* Broadcast cell data to all other processes */
        MPI_Bcast(&Ncells, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(me != 0) {
            cells = (Cell*) malloc(Ncells*sizeof(Cell));
            if(cells == NULL) {
                fprintf(stderr, "ReadCells: could not allocate memory for cells on process %d\n", me);
                opsec_abort(1);
            }
        }
        MPI_Bcast(cells, Ncells, cell_datatype, 0, MPI_COMM_WORLD);

        MPI_Type_free(&cell_datatype);
    }
#endif

    if(Ncells_) *Ncells_ = Ncells;
    return cells;
}
