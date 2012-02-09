#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <mpi.h>

#include "Cell.h"
#include "Matrix.h"
#include "abn.h"
#include "cfg.h"

void print_usage(FILE* f, char* progname) {
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0)
        fprintf(f,
"Usage: %s [SWITCHES] [name=value ...]\n"
"Switches:\n"
"  -h            Display this help\n"
"  -c FILE       Read additional configuration options from FILE\n"
"Configuration options:\n"
"  cellfile=FILE Read list of cells from FILE (required)\n"
"  specfile=FILE Read KL spectrum from FILE (required)\n",
                progname);
}

void report_memory_usage() {
    printf("Memory usage from /proc/%d/statm: ", getpid());
    fflush(stdout);
    char cmd[256];
    snprintf(cmd, sizeof(cmd), "cat /proc/%d/statm", getpid());
    system(cmd);
}

Cell* ReadCells(MPI_Comm comm, const char* cellfile, int &Ncells) {
    int nprocs, myid, root;
    Cell* cells;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    root = 0;

    /* Have root process read in cells */
    if(myid == root) {
        size_t n, size;
        char endian;
        char fmt[ABN_MAX_FORMAT_LENGTH];
        Config cellopts = cfg_new();

        /* Read basis cells from file */
        FILE* fcells = fopen(cellfile, "r");
        if(fcells == NULL) {
            fprintf(stderr, "ReadCells: could not read cells from '%s'\n", cellfile);
            MPI_Abort(comm, 1);
        }
        if(abn_read(fcells, (void**) &cells, &n, &size, &endian, fmt, cellopts) != 0) {
            report_memory_usage();
            fprintf(stderr, "ReadCells: error reading cells from '%s'\n", cellfile);
            MPI_Abort(comm, 1);
        }
        fclose(fcells);

        Ncells = cfg_get_int(cellopts, "Ncells");
        if((int)n != Ncells)
            fprintf(stderr, "Warning: consistency check fail for Ncells: %d != %d\n", Ncells, (int)n);
        if(size != sizeof(Cell))
            fprintf(stderr, "Warning: consistency check fail for sizeof(Cell): %d != %d\n", (int)sizeof(Cell), (int)size);

        /* Debugging... */
        printf("Read %d basis cells from '%s'\n", Ncells, cellfile);
    }
    MPI_Bcast(&Ncells, 1, MPI_INT, root, comm);

    /* Define MPI datatype for cells (this must be updated if struct Cell changes!) */
    MPI_Datatype cell_datatype;
    int blocklengths[2] = { 2, 8 };
    MPI_Aint displacements[2] = { 0, 8 };
    MPI_Datatype types[2] = { MPI_INT, MPI_DOUBLE };
    MPI_Type_create_struct(2, blocklengths, displacements, types, &cell_datatype);
    MPI_Type_commit(&cell_datatype);

    /* Broadcast cell data to all other processes */
    if(myid != root) {
        cells = (Cell*) malloc(Ncells*sizeof(Cell));
        if(cells == NULL) {
            fprintf(stderr, "ReadCells: could not allocate memory for cells on process %d\n", myid);
            MPI_Abort(comm, 1);
        }
    }
    MPI_Bcast(cells, Ncells, cell_datatype, root, comm);
    return cells;
}

double* ReadSpectrum(const char* specfile, int* n_, int* nev_) {
    int n, nev;
    double* modes = NULL;
    Config opts = cfg_new();

    FILE* fspec = fopen(specfile, "r");
    if(fspec == NULL) {
        fprintf(stderr, "ReadSpectrum: could not open file '%s'\n", specfile);
        return NULL;
    }

    int err = abn_read(fspec, (void**) &modes, NULL, NULL, NULL, NULL, opts);
    fclose(fspec);

    if(!err && !cfg_has_keys(opts, "n,nev", ",")) {
        fprintf(stderr, "ReadSpectrum: missing options\n");
        err = 1;
    }
    else {
        n = cfg_get_int(opts, "n");
        nev = cfg_get_int(opts, "nev");
    }
    cfg_destroy(opts);

    if(err) {
        fprintf(stderr, "ReadSpectrum: error reading from '%s'\n", specfile);
        perror("system error");
        free(modes);
        return NULL;
    }
    else {
        if(n_) *n_ = n;
        if(nev_) *nev_ = nev;
        return modes;
    }
}

int ReadSpectrum(MPI_Comm comm, const char* specfile, int& Nmodes, int& Ncells, Matrix& B, Matrix& C) {
    int nprocs, myid, root;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    root = 0;

    /* Read in spectrum from file on the root process */
    if(myid == root) {
        int n, nev;
        double* modes = ReadSpectrum(specfile, &n, &nev);
        if(modes == NULL) {
            Ncells = -1;
            Nmodes = -1;
        }
        else {
            Ncells = n;
            Nmodes = nev;
            B.resize(Nmodes, Ncells);
            C.resize(Nmodes, Nmodes);

            /* Eigenvalues of the cell-basis signal matrix become the diagonal elements
             * of the covariance matrix in the KL basis.  Eigenvectors become the rows
             * of the projection matrix. */
            for(int i = 0; i < Nmodes; i++) {
                for(int j = 0; j < Nmodes; j++)
                    C(i,j) = (i == j) ? 1 + modes[i*(Ncells+1)] : 0;
                for(int a = 0; a < Ncells; a++)
                    B(i,a) = modes[i*(Ncells+1) + a + 1];
            }
            free(modes);
        }
    }

    /* Broadcast spectrum to other processes */
    MPI_Bcast(&Nmodes, 1, MPI_INT, root, comm);
    MPI_Bcast(&Ncells, 1, MPI_INT, root, comm);
    if(Nmodes > 0) {
        if(myid != root) {
            B.resize(Nmodes, Ncells);
            C.resize(Nmodes, Nmodes);
        }
        MPI_Bcast(B.data, Nmodes*Ncells, MPI_DOUBLE, root, comm);
        MPI_Bcast(C.data, Nmodes*Nmodes, MPI_DOUBLE, root, comm);
        return 0;
    }
    else {
        fprintf(stderr, "ReadSpectrum: spectrum not read\n");
        return 1;
    }
}


int main(int argc, char* argv[]) {
    int nprocs, myid, root;
    int Ncells, Nmodes;
    Cell* cells;

    Config cfg = cfg_new();

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myid);
    root = nprocs - 1;

    /* Parse command line switches */
    int opt;
    const char* optstring = "hc:";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            print_usage(stdout, argv[0]);
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            print_usage(stderr, argv[0]);
            MPI_Abort(comm, 1);
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Debugging... */
    if(myid == root) {
        printf("# Config options\n");
        cfg_write(cfg, stdout);
    }

    if(myid == 0) {
        printf("Testing malloc...\n");
        char* test = (char*) malloc(72000*sizeof(char));
        if(test == NULL) {
            fprintf(stderr, "Could not allocate a measly 72000 bytes!\n");
            perror("system error:");
        }
        free(test);
    }

    /* Make sure all the necessary options are provided */
    if(!cfg_has_keys(cfg, "cellfile,specfile", ",")) {
        print_usage(stderr, argv[0]);
        MPI_Abort(comm, 1);
    }
    const char* cellfile = cfg_get(cfg, "cellfile");
    const char* specfile = cfg_get(cfg, "specfile");

    /* Read cells from file */
    if(myid == 0) printf("Reading cells from '%s' ...\n", cellfile);
    cells = ReadCells(comm, cellfile, Ncells);
    if(cells == NULL) {
        if(myid == 0) fprintf(stderr, "Failed to read cells.\n");
        if(myid == 0) perror("system error");
        MPI_Abort(comm, 1);
    }
    else {
        if(myid == 0) printf("Read %d cells.\n", Ncells);
        if(myid == 0) printf("  cells[0] = { %d, %d, %g, %g, %g, %g, %g, %g, %g, %g }\n", cells[0].a, cells[0].G, cells[0].rmin, cells[0].rmax, cells[0].mumin, cells[0].mumax, cells[0].phimin, cells[0].phimax, cells[0].Veff, cells[0].Nbar);
    }

    /* Read KL mode coefficients */
    Matrix B;   // cell -> KL mode projection matrix
    Matrix C;   // full covariance matrix in KL basis (diagonal)
    if(myid == 0) printf("Reading KL spectrum from '%s' ...\n", specfile);
    if(ReadSpectrum(comm, specfile, Nmodes, Ncells, B, C) != 0) {
        if(myid == 0) fprintf(stderr, "Failed to read spectrum.\n");
        if(myid == 0) perror("system error");
        MPI_Abort(comm, 1);
    }
    else {
        if(myid == 0) printf("Read spectrum: Nmodes = %d, Ncells = %d\n", Nmodes, Ncells);
    }

    /* Clean up nicely */
    free(cells);
    MPI_Finalize();

    return 0;
}
