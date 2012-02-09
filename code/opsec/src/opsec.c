#ifdef HAVE_CONFIG_H
#  include <opsec_config.h>
#endif

#include <getopt.h>
#ifdef OPSEC_USE_MPI
#  include <mpi.h>
#endif
#ifdef _OPENMP
#  include <omp.h>
#endif

#include "opsec.h"

/* OPSEC globals */
int opsec_verbose = 0;


Config opsec_init(int argc, char * const argv[], const char *usage) {
    Config cfg = cfg_new();
    int i;

    int opt;
    const char *optstring = "hc:";

    /* Parse command line switches */
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            fputs(usage, stdout);
            opsec_exit(0);
            break;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fputs(usage, stdout);
            opsec_exit(1);
            break;
        }
    }

    /* Parse additional command line options */
    for(i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Process special values */
#ifdef _OPENMP
    if(cfg_has_key(cfg, "num_threads")) {
        int num_threads = cfg_get_int(cfg, "num_threads");
        omp_set_num_threads(num_threads);
    }
#endif

    return cfg;
}

void opsec_exit(int status) {
    /* This routine should work for both MPI and non-MPI programs, including
     * those that are compiled with OPSEC_USE_MPI present, but which never call
     * MPI_Init().  So we have to be a bit careful. */
    int me = 0;

#ifdef OPSEC_USE_MPI
    /* Decide if MPI_Init() was ever called. */
    int mpi;
    MPI_Initialized(&mpi);

    /* Decide if we're the root process. */
    if(mpi)
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    if(me == 0) {
        printf("OPSEC is now exiting.\n\n");
        fflush(stdout);
        fflush(stderr);
    }

#ifdef OPSEC_USE_MPI
    if(mpi)
        MPI_Finalize();
#endif

    exit(status);
}

void opsec_abort(int status) {
    /* This routine should work for both MPI and non-MPI programs, including
     * those that are compiled with OPSEC_USE_MPI present, but which never call
     * MPI_Init().  So we have to be a bit careful. */
    int me = 0;

#ifdef OPSEC_USE_MPI
    /* Decide if MPI_Init() was ever called. */
    int mpi;
    MPI_Initialized(&mpi);

    /* Decide if we're the root process. */
    if(mpi)
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    if(me == 0) {
        fprintf(stderr, "OPSEC is aborting!!!.\n\n");
        fflush(stdout);
        fflush(stderr);
    }

#ifdef OPSEC_USE_MPI
    if(mpi)
        MPI_Abort(MPI_COMM_WORLD, status);
#endif

    abort();
}
