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

int opsec_info(const char* format, ...) {
    int nchars = 0;
    va_list ap;
    va_start(ap, format);
    nchars = opsec_vfprintf(stdout, format, ap);
    va_end(ap);
    return nchars;
}

int opsec_debug(const char* format, ...) {
    int nchars = 0;
#ifdef OPSEC_DEBUG
    va_list ap;
    va_start(ap, format);
    nchars = opsec_vfprintf(stdout, format, ap);
    va_end(ap);
#endif
    return nchars;
}

int opsec_warn(const char* format, ...) {
    int nchars = 0;
#ifdef OPSEC_WARN
    va_list ap;
    va_start(ap, format);
    nchars = opsec_vfprintf(stderr, format, ap);
    va_end(ap);
#endif
    return nchars;
}

int opsec_error(const char* format, ...) {
    int nchars = 0;
    va_list ap;
    va_start(ap, format);
    nchars = opsec_vfprintf(stderr, format, ap);
    va_end(ap);
    return nchars;
}

int opsec_vfprintf(FILE* stream, const char *format, va_list ap) {
    int nchars = 0;

    /* This routine should work for both MPI and non-MPI programs, including
     * those that are compiled with OPSEC_USE_MPI present, but which never call
     * MPI_Init().  So we have to be a bit careful. */
    int nprocs = 1, me = 0;

#ifdef OPSEC_USE_MPI
    /* Decide if MPI_Init() was ever called. */
    int mpi;
    MPI_Initialized(&mpi);

    /* Decide if we're the root process. */
    if(mpi) {
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
    }
#endif
    if(nprocs > 1)
        nchars = fprintf(stream, "[proc %d] ", me);
    nchars += vfprintf(stream, format, ap);
    fflush(stream);
    return nchars;
}

void* opsec_debug_malloc(size_t size, const char* file, int line) {
    void* ptr = malloc(size);
    if(!ptr) {
        opsec_error("debug_malloc: failed to allocate %zd bytes at line %d of %s\n", size, line, file);
        opsec_abort(1);
    }
    return ptr;
}
void* opsec_debug_calloc(size_t nmemb, size_t size, const char* file, int line) {
    void* ptr = calloc(nmemb, size);
    if(!ptr) {
        opsec_error("debug_calloc: failed to allocate %zd bytes at line %d of %s\n", nmemb*size, line, file);
        opsec_abort(1);
    }
    return ptr;
}
