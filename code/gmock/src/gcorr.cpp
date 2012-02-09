/* gcorr
 *
 * Generate a Gaussian random density field and it's associated linear velocity
 * field, and compute the real- and redshift-space correlation functions from
 * it.
 *
 * Algorithm steps:
 * 1. Generate density and velocity fields
 * 2. Sample particles
 * 3. Compute correlation function
 * 4. Repeat 2 and 3 until correlation function converges
 *
 * Notation:
 * - nx,ny,nz are the number of grid cells in each direction within the
 *   computational box
 * - nxloc is the number of grid cells in the x direction within the current
 *   slab
 * - ix,iy,iz index a particular grid cell within the computational box
 * - ixloc indexes grid cells within the current slab
 *
 * - mx,my,mz are the number of chaining mesh cells in each direction within
 *   the computational box
 * - mxloc is the number of chaining mesh cells in the x direction within the
 *   current slab
 * - jx,jy,jz index a particular chaining mesh cell within the computational
 *   box
 * - jxloc indexes chaining mesh cells within the current slab
 *
 * - hx,hy,hz index grid cells within the current chaining mesh cell
 *
 *
 * TODO:
 * - remove restrictions on N, M, and nprocs
 */

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include "abn.h"
#include "cfg.h"
#include "estimator.h"
#include "grid.h"
#include "particle.h"
#include "rng.h"
#include "sampling.h"
#include "Spline.h"

const unsigned int DEFAULT_SEED = 1776;

/* Information describing a computational slab (i.e. the portion of the
 * computational box that belongs to a single process) */
struct Slab {
    double L;                   // dimensions of computational box
    int N;                      // number of grid cells in computational box (in each direction)
    int M;                      // number of chaining mesh cells in computational box (in each direction)
    int nxloc, ny, nz, ixmin;   // number of grid cells in slab, and position of slab within grid
    int mxloc, my, mz, jxmin;   // number of mesh cells in slab, and position of slab within chaining mesh
    Grid rho, vx, vy, vz;       // density and velocity fields
};

/* Inverse CIC: return the value of the field at the point (x,y,z) \in [0,1]^3 */
static inline double icic(double a000, double a001, double a010, double a011,
                          double a100, double a101, double a110, double a111,
                          double x, double y, double z)
{
    return a000*(1-x)*(1-y)*(1-z) + a001*(1-x)*(1-y)*z + a010*(1-x)*y*(1-z) + a011*(1-x)*y*z
         + a100*x*(1-y)*(1-z)     + a101*x*(1-y)*z     + a110*x*y*(1-z)     + a111*x*y*z;
}

/* Populate cell (jxloc,jy,jz) within the chaining mesh with particles, drawn
 * from a number density field rho. */
void populate_cell(int jxloc, int jy, int jz, ChainingMesh& mesh, Slab& slab) {
    double L = slab.L;
    int N = slab.N;
    int M = slab.M;
    int nxloc = slab.nxloc, ny = slab.ny, nz = slab.nz, ixmin = slab.ixmin;
    int mxloc = slab.mxloc, my = slab.my, mz = slab.mz, jxmin = slab.jxmin;
    assert((nxloc % mxloc) == 0 && (ny % my) == 0 && (nz % mz) == 0);
    double dL = L/N;            // grid cell size
    double dV = dL*dL*dL;       // grid cell volume
    Grid& rho = slab.rho;
    Grid& vx = slab.vx;
    Grid& vy = slab.vy;
    Grid& vz = slab.vz;

    int jx = jxmin + jxloc;
    Cell& cell = mesh.cells[(jxloc*my + jy)*mz + jz];
    cell.xmin = jx*L/M;
    cell.xmax = (jx+1)*L/M;
    cell.ymin = jy*L/M;
    cell.ymax = (jy+1)*L/M;
    cell.zmin = jz*L/M;
    cell.zmax = (jz+1)*L/M;
    cell.start = (int)mesh.particles.size();
    cell.count = 0;

    /* Iterate over grid cells within this chaining mesh cell */
    for(int hx = 0; hx < nxloc/mxloc; hx++) {
        int ixloc = hx + jxloc*(nxloc/mxloc);
        int ix = ixloc + ixmin;
        for(int hy = 0; hy < ny/my; hy++) {
            int iy = hy + jy*(ny/my);
            for(int hz = 0; hz < nz/mz; hz++) {
                int iz = hz + jz*(nz/mz);

                /* Compute mean density over the grid cell */
                double a000 = rho(ixloc,iy,iz);
                double a001 = rho(ixloc,iy,(iz+1)%nz);
                double a010 = rho(ixloc,(iy+1)%ny,iz);
                double a011 = rho(ixloc,(iy+1)%ny,(iz+1)%nz);
                double a100 = rho(ixloc+1,iy,iz);
                double a101 = rho(ixloc+1,iy,(iz+1)%nz);
                double a110 = rho(ixloc+1,(iy+1)%ny,iz);
                double a111 = rho(ixloc+1,(iy+1)%ny,(iz+1)%nz);
                double rhobar = (a000 + a001 + a010 + a011 + a100 + a101 + a110 + a111)/8;

                /* Decide how many new particles to add to this volume element */
                int numnew = (int)rng_poisson(rhobar*dV);
                double fx, fy, fz;      // position of particle within grid cell, (fx,fy,fz) \in [0,1]^3
                Particle p;
                for(int i = 0; i < numnew; i++) {
                    /* Choose where within the grid cell to place the particle */
                    sample3d(a000, a001, a010, a011, a100, a101, a110, a111, fx, fy, fz);
                    p.x = (ix + fx)*dL;
                    p.y = (iy + fy)*dL;
                    p.z = (iz + fz)*dL;

                    /* Compute the velocity at that position */
                    p.vx = icic(vx(ixloc,iy,iz), vx(ixloc,iy,(iz+1)%N), vx(ixloc,(iy+1)%N,iz), vx(ixloc,(iy+1)%N,(iz+1)%N), vx(ixloc+1,iy,iz), vx(ixloc+1,iy,(iz+1)%N), vx(ixloc+1,(iy+1)%N,iz), vx(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);
                    p.vy = icic(vy(ixloc,iy,iz), vy(ixloc,iy,(iz+1)%N), vy(ixloc,(iy+1)%N,iz), vy(ixloc,(iy+1)%N,(iz+1)%N), vy(ixloc+1,iy,iz), vy(ixloc+1,iy,(iz+1)%N), vy(ixloc+1,(iy+1)%N,iz), vy(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);
                    p.vz = icic(vz(ixloc,iy,iz), vz(ixloc,iy,(iz+1)%N), vz(ixloc,(iy+1)%N,iz), vz(ixloc,(iy+1)%N,(iz+1)%N), vz(ixloc+1,iy,iz), vz(ixloc+1,iy,(iz+1)%N), vz(ixloc+1,(iy+1)%N,iz), vz(ixloc+1,(iy+1)%N,(iz+1)%N), fx, fy, fz);

                    /* Add the particle */
                    mesh.particles.push_back(p);
                    cell.count += 1;
                }
            }
        }
    }
}

/* Sample data particles from within this slab, with number density field given
 * by rho and velocity field given by (vx,vy,vz). */
void populate_mesh_data(ChainingMesh& mesh, Slab& slab) {
    int mxloc = slab.mxloc, my = slab.my, mz = slab.mz;
    mesh.cells.resize(mxloc*my*mz);
    mesh.particles.clear();

    /* Iterate over chaining mesh cells within slab */
    for(int jxloc = 0; jxloc < mxloc; jxloc++)
        for(int jy = 0; jy < my; jy++)
            for(int jz = 0; jz < mz; jz++)
                populate_cell(jxloc, jy, jz, mesh, slab);
}

/* Sample random particles from within this slab, with uniform number density
 * rhobar. */
void populate_mesh_random(ChainingMesh& mesh, Slab& slab, double rhobar) {
    double L = slab.L;
    int M = slab.M;
    int mxloc = slab.mxloc, my = slab.my, mz = slab.mz, jxmin = slab.jxmin;
    double dL = L/M;            // size of individual chaining mesh cell
    double dV = dL*dL*dL;       // volume of individual chaining mesh cell

    mesh.cells.resize(mxloc*my*mz);
    mesh.particles.clear();

    /* Iterate over chaining mesh cells within slab */
    for(int jxloc = 0; jxloc < mxloc; jxloc++) {
        int jx = jxmin + jxloc;
        for(int jy = 0; jy < my; jy++) {
            for(int jz = 0; jz < mz; jz++) {
                Cell& cell = mesh.cells[(jxloc*my + jy)*mz + jz];
                cell.xmin = jx*L/M;
                cell.xmax = (jx+1)*L/M;
                cell.ymin = jy*L/M;
                cell.ymax = (jy+1)*L/M;
                cell.zmin = jz*L/M;
                cell.zmax = (jz+1)*L/M;
                cell.count = 0;
                cell.start = (int)mesh.particles.size();

                /* Decide how many new particles to add to this volume element */
                int numnew = (int)rng_poisson(rhobar*dV);
                Particle p;
                for(int i = 0; i < numnew; i++) {
                    /* Choose where within the chaining mesh cell to place the particle */
                    p.x = (jx + rng_uniform())*dL;
                    p.y = (jy + rng_uniform())*dL;
                    p.z = (jz + rng_uniform())*dL;

                    p.vx = p.vy = p.vz = 0;

                    /* Add the particle */
                    mesh.particles.push_back(p);
                    cell.count++;
                }
            }
        }
    }
}

static bool datatype_initialized = false;
static MPI_Datatype mpi_cell_datatype;

void send_mesh(MPI_Comm comm, ChainingMesh& mesh, int dest) {
    if(!datatype_initialized) {
        /* Establish MPI datatype for Cell */
        int blocklengths[2] = { 6, 2 };
        MPI_Aint displacements[2] = { 0, 6*8 };
        MPI_Datatype types[2] = { MPI_DOUBLE, MPI_INT };
        MPI_Type_create_struct(2, blocklengths, displacements, types, &mpi_cell_datatype);
        MPI_Type_commit(&mpi_cell_datatype);
        datatype_initialized = true;
    }

    int ncells = (int)mesh.cells.size();
    int nparticles = (int)mesh.particles.size();
    MPI_Send(&ncells, 1, MPI_INT, dest, 0, comm);
    MPI_Send(&nparticles, 1, MPI_INT, dest, 0, comm);
    MPI_Send(&mesh.cells[0], ncells, mpi_cell_datatype, dest, 0, comm);
    MPI_Send(&mesh.particles[0], 6*nparticles, MPI_FLOAT, dest, 0, comm);
}

void recv_mesh(MPI_Comm comm, ChainingMesh& mesh, int source) {
    if(!datatype_initialized) {
        /* Establish MPI datatype for Cell */
        int blocklengths[2] = { 6, 2 };
        MPI_Aint displacements[2] = { 0, 6*8 };
        MPI_Datatype types[2] = { MPI_DOUBLE, MPI_INT };
        MPI_Type_create_struct(2, blocklengths, displacements, types, &mpi_cell_datatype);
        MPI_Type_commit(&mpi_cell_datatype);
        datatype_initialized = true;
    }

    int ncells, nparticles;
    MPI_Recv(&ncells, 1, MPI_INT, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&nparticles, 1, MPI_INT, source, 0, comm, MPI_STATUS_IGNORE);
    mesh.cells.resize(ncells);
    mesh.particles.resize(nparticles);
    MPI_Recv(&mesh.cells[0], ncells, mpi_cell_datatype, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&mesh.particles[0], 6*nparticles, MPI_FLOAT, source, 0, comm, MPI_STATUS_IGNORE);
}

void sample_particles(MPI_Comm comm, Binner& binner, double smax, Slab& slab, double rhobar, double& ndata, double& nrandom) {
    int nprocs, me;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    ChainingMesh pdata, prandom, other;

    /* Fill chaining meshes with data and random particles */
    populate_mesh_data(pdata, slab);
    populate_mesh_random(prandom, slab, rhobar);

    /* Find maximum velocity of data particles */
    double vmax = 0.;
    for(vector<Particle>::const_iterator p = pdata.particles.begin(); p != pdata.particles.end(); p++)
        vmax = fmax(vmax, sqrt(p->vx*p->vx + p->vy*p->vy + p->vz*p->vz));
    MPI_Allreduce(MPI_IN_PLACE, &vmax, 1, MPI_DOUBLE, MPI_MAX, comm);

//        printf("  vmax = %g\n", vmax);
//        printf("  ndata = %u, nrandom = %u\n", pdata.particles.size(), prandom.particles.size());

    /* Count data-data pairs */
    binner.self_correlate_dd(pdata, smax + 2*vmax);
    for(int hop = 1; hop <= nprocs/2; hop++) {
        /* Send our data particles to another process, and receive another
         * process' data particles in kind */
        for(int source = 0; source < nprocs; source++) {
            int dest = (source + hop) % nprocs;
            if(me == source)
                send_mesh(comm, pdata, dest);
            else if(me == dest)
                recv_mesh(comm, other, source);
        }

        /* Now correlate our data with theirs (on the last hop, if nprocs is
         * even, only half the processes need to do work) */
        if(2*hop != nprocs || me < nprocs/2)
            binner.correlate_dd(pdata, other, smax + 2*vmax);

        MPI_Barrier(comm);
    }

    /* Count random-random pairs */
    binner.self_correlate_rr(prandom, smax);
    for(int hop = 1; hop <= nprocs/2; hop++) {
        /* Send our random particles to another process, and receive another
         * process' random particles in kind */
        for(int source = 0; source < nprocs; source++) {
            int dest = (source + hop) % nprocs;
            if(me == source)
                send_mesh(comm, prandom, dest);
            else if(me == dest)
                recv_mesh(comm, other, source);
        }

        /* Now correlate our randoms with theirs (on the last hop, if nprocs is
         * even, only half the processes need to do work) */
        if(2*hop != nprocs || me < nprocs/2)
            binner.correlate_rr(prandom, other, smax);
    }

    /* Count data-random pairs */
    binner.correlate_dr(pdata, prandom, smax + vmax);
    for(int hop = 1; hop < nprocs; hop++) {
        /* Send our data particles to another process, and receive another
         * process' data particles in kind (we transfer the data and not the
         * randoms because there are typically more of the latter) */
        for(int source = 0; source < nprocs; source++) {
            int dest = (source + hop) % nprocs;
            if(me == source)
                send_mesh(comm, pdata, dest);
            else if(me == dest)
                recv_mesh(comm, other, source);
        }

        /* Now correlate our randoms with their data */
        binner.correlate_dr(other, prandom, smax + vmax);
    }

    /* Determine the number of data and random particles sampled */
    ndata = (double)pdata.particles.size();
    nrandom = (double)prandom.particles.size();
    MPI_Allreduce(MPI_IN_PLACE, &ndata, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &nrandom, 1, MPI_DOUBLE, MPI_SUM, comm);
}

/* Check whether arrays a and b are similar (for each element, either the
 * absolute difference is less than epsabs, or the relative difference is
 * less than epsrel) */
bool similar(const array& a, const array& b, double epsrel, double epsabs) {
    assert(epsrel > 0 && epsabs > 0);
    int N = (int)a.size();
    double d, t;
    for(int i = 0; i < N; i++) {
        d = fabs(b[i] - a[i]);
        t = fabs(a[i]) + epsabs;      // + epsabs is just to prevent divide-by-zero
        if(d > epsabs && d/t > epsrel)
            return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    int nprocs = 1, me = 0;     // MPI parameters (default is compatible with non-MPI build)

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    /* Parse command line and configuration options */
    Config cfg = cfg_new();
    int opt;
    while((opt = getopt(argc, argv, "hc:")) != -1) {
        switch(opt) {
        case 'h':
            printf("Usage: gcorr -c cfgfile\n");
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            fprintf(stderr, "Invalid command line options.  Try 'gcorr -h'.\n");
            return 1;
        }
    }
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Debugging */
    if(me == 0) {
        printf("Config options:\n");
        cfg_write(cfg, stdout);
        printf("\n");
    }

//    if(!cfg_has_keys(cfg, "N,L,f,b,sigma,rhobar,outfile", ",")) {
    if(!cfg_has_keys(cfg, "N,L,f,b,Rsmooth,rhobar,outfile", ",")) {
        fprintf(stderr, "Missing configuration options.\n");
        return 1;
    }
    const char* pkfile = cfg_get(cfg, "pkfile");
    int N = cfg_get_int(cfg, "N");
    int M = cfg_has_key(cfg, "M") ? cfg_get_int(cfg, "M") : nprocs;
    double L = cfg_get_double(cfg, "L");
    double f = cfg_get_double(cfg, "f");
    double b = cfg_get_double(cfg, "b");
    double rhobar = cfg_get_double(cfg, "rhobar");
//    double sigma = cfg_get_double(cfg, "sigma");
    double Rsmooth = cfg_get_double(cfg, "Rsmooth");
    const char* outfile = cfg_get(cfg, "outfile");
    int maxsamp = cfg_has_key(cfg, "maxsamp") ? cfg_get_int(cfg, "maxsamp") : 1000;
    double epsrel0 = cfg_has_key(cfg, "epsrel0") ? cfg_get_double(cfg, "epsrel0") : 0.01;
    double epsabs0 = cfg_has_key(cfg, "epsabs0") ? cfg_get_double(cfg, "epsabs0") : 1e-5;
    double epsrel2 = cfg_has_key(cfg, "epsrel2") ? cfg_get_double(cfg, "epsrel2") : 0.1;
    double epsabs2 = cfg_has_key(cfg, "epsabs2") ? cfg_get_double(cfg, "epsabs2") : 1e-5;
    double epsrel4 = cfg_has_key(cfg, "epsrel4") ? cfg_get_double(cfg, "epsrel4") : 1.;
    double epsabs4 = cfg_has_key(cfg, "epsabs4") ? cfg_get_double(cfg, "epsabs4") : 1e-5;
    int nsbins = cfg_has_key(cfg, "nsbins") ? cfg_get_int(cfg, "nsbins") : 39;
    int nmubins = cfg_has_key(cfg, "nmubins") ? cfg_get_int(cfg, "nmubins") : 64;
    double smin = cfg_has_key(cfg, "smin") ? cfg_get_double(cfg, "smin") : 3.;
    double smax = cfg_has_key(cfg, "smax") ? cfg_get_double(cfg, "smax") : 120.;

    /* Impose parameter restrictions: things are way easier this way */
    if((N % M) != 0 || (M % nprocs) != 0) {
        fprintf(stderr, "Error: M must be a multiple of nprocs, and N must be a multiple of M (N = %d, M = %d, nprocs = %d)\n", N, M, nprocs);
        return 1;
    }

    /* Initialize slab information */
    Slab slab;
    slab.N = N;
    slab.M = M;
    slab.L = L;
    slab.nxloc = N/nprocs;
    slab.ny = slab.nz = N;
    slab.ixmin = me*slab.nxloc;
    slab.mxloc = M/nprocs;
    slab.my = slab.mz = M;
    slab.jxmin = me*slab.mxloc;
    Grid& rho = slab.rho;
    Grid& vx = slab.vx;
    Grid& vy = slab.vy;
    Grid& vz = slab.vz;

    /* Predict how much memory per node this is going to take, and halt if we
     * don't think there's going to be enough */
    // TODO

    /* Make sure we can open the output file */
    FILE* fout = NULL;
    if(me == 0) {
        fout = fopen(outfile, "w");
        if(!fout) {
            fprintf(stderr, "gcorr: cannot open '%s' for writing.\n", outfile);
            return 1;
        }
    }

    /* Seed random number generator.  All processes choose the same seed: this
     * is important for generating the density field. */
    unsigned long seed = cfg_has_key(cfg, "seed") ? cfg_get_uint(cfg, "seed") : DEFAULT_SEED;
    rng_init(seed);

    if(me == 0)
        printf("Loading power spectrum\n");

    /* Load power spectrum from file */
//    Spline P = LinearSpline(pkfile);

    /* Model power spectrum by P(k) \propto k e^{-kR} */
    const int Nk = 1024;
    const double kmin = 0., kmax = 1.;
    const double R = 50.;
    array kk(Nk), pk(Nk);
    for(int i = 0; i < Nk; i++) {
        kk[i] = kmin + i*(kmax - kmin)/(Nk-1);
        pk[i] = kk[i] * exp(-kk[i]*R);
    }
    Spline P = LinearSpline(kk, pk);

    if(me == 0)
        printf("Preparing grid\n");

    /* Allocate grids */
    rho.initialize(comm, N, N, N);
    vx.initialize(comm, N, N, N);
    vy.initialize(comm, N, N, N);
    vz.initialize(comm, N, N, N);

    /* Make sure the FFTW memory allocation is as expected */
    assert(slab.nxloc == rho.nxloc && slab.ixmin == rho.ixmin);
    int nxloc = slab.nxloc;

    /*! Note that initially the Grid rho actually represents the matter density
     *! contrast (delta).  After rescaling below, it represents the galaxy
     *! number density normalized to have a mean of rhobar. */

    if(me == 0)
        printf("Generating density and velocity fields\n");

    generate_fields(N, L, f, P, Rsmooth, rho, vx, vy, vz);

#if 0
    /* Rescale density contrast to have desired sigma */
    {
        double var = 0.;
        for(int ixloc = 0; ixloc < nxloc; ixloc++)
            for(int iy = 0; iy < N; iy++)
                for(int iz = 0; iz < N; iz++)
                    var += pow(rho(ixloc,iy,iz), 2);

        if(me == 0)
            MPI_Reduce(MPI_IN_PLACE, &var, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        else
            MPI_Reduce(&var, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

        double norm = sigma*sqrt(N*N*N/var);
        MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, comm);
        if(me == 0)
            printf("  norm = %g\n", norm);

        for(int ixloc = 0; ixloc < nxloc; ixloc++) {
            for(int iy = 0; iy < N; iy++) {
                for(int iz = 0; iz < N; iz++) {
                    rho(ixloc,iy,iz) *= norm;
                    vx(ixloc,iy,iz) *= norm;
                    vy(ixloc,iy,iz) *= norm;
                    vz(ixloc,iy,iz) *= norm;
                }
            }
        }
    }
#endif

    /* Convert density contrast to proper number density */
    double numneg = 0;
    for(int ixloc = 0; ixloc < nxloc; ixloc++) {
        for(int iy = 0; iy < N; iy++) {
            for(int iz = 0; iz < N; iz++) {
                if(rho(ixloc,iy,iz) < -1/b)
                    numneg += 1;
                rho(ixloc,iy,iz) = rhobar * fmax(0, 1 + b*rho(ixloc,iy,iz));
            }
        }
    }
    printf("Process %d: density field %g percent negative\n", me, 100*numneg/(nxloc*N*N));

    if(me == 0)
        printf("Transferring boundary layers\n");

    /* Use the extra memory allocated by Grid (nxloc+1 instead of nxloc) to
     * store the grid values from the adjacent process, process me+1.  At the
     * same time, we need to send our boundary grid values to process me-1. */
    int dest = (nprocs+me-1) % nprocs;
    int source = (me+1) % nprocs;
    int count = N*2*(N/2+1);
    assert(sizeof(fftw_real) == sizeof(double));
    MPI_Sendrecv(&rho(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &rho(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vx(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vx(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vy(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vy(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&vz(0,0,0), count, MPI_DOUBLE, dest, 0,
                 &vz(nxloc,0,0), count, MPI_DOUBLE, source, 0, comm, MPI_STATUS_IGNORE);

    /* Now that the grid has been generated, reseed the random number generator
     * with a different seed on each process.  This isn't too important, as the
     * different density field on each slab will cause the random numbers to
     * diverge, but it's best to be safe. */
    rng_init(seed + 97*me);     // 97 is arbitrary

    PlaneParallelBinner binner(L, nsbins, nmubins, smin, smax);
    PlaneParallelEstimator estimator(nsbins, nmubins);

    /* Current and previous estimates for Legendre moments \xi_\ell */
    array xi0, xi2, xi4, xi0_last, xi2_last, xi4_last;
    xi0_last = xi2_last = xi4_last = array::zeros(nsbins);

    int nconv = 0;       // number of consecutive samplings for which convergence tests pass
    for(int isamp = 0; nconv < 5 && isamp < maxsamp; isamp++) {
        if(me == 0)
            printf("Starting sampling number %d (max %d)\n", isamp+1, maxsamp);

        double ndata, nrandom;  // number of new particles added for this sampling
        sample_particles(comm, binner, smax, slab, rhobar, ndata, nrandom);

        /* Update \xi(s,\mu) estimate, and check for convergence */
        if(me == 0) {
            MPI_Reduce(MPI_IN_PLACE, &ndata, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(MPI_IN_PLACE, &nrandom, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.dd[0], &estimator.dd[0], nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.dr[0], &estimator.dr[0], nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.rr[0], &estimator.rr[0], nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);

            /* Update total pair counts */
            estimator.ndd += ndata*(ndata+1.)/2.;
            estimator.ndr += ndata*nrandom;
            estimator.nrr += nrandom*(nrandom+1.)/2.;

            /* Update Legendre moment estimates */
            xi0 = estimator.legendre(0);
            xi2 = estimator.legendre(2);
            xi4 = estimator.legendre(4);

            /* Check for convergence */
            bool c = similar(xi0, xi0_last, epsrel0, epsabs0)
                  && similar(xi2, xi2_last, epsrel2, epsabs2)
                  && similar(xi4, xi4_last, epsrel4, epsabs4);
            if(c)
                nconv++;
            else
                nconv = 0;

//            double dabs = max(abs(xi0-xi0_last));
//            double drel = max(abs(xi0-xi0_last)/(abs(xi0)+epsabs));
//            printf("  xi_0: dabs = %g, drel = %g\n", dabs, drel);
//            converged = (dabs < epsabs) || (drel < epsrel);
//            dabs = max(abs(xi2-xi2_last));
//            drel = max(abs(xi2-xi2_last)/(abs(xi2)+epsabs));
//            printf("  xi_2: dabs = %g, drel = %g\n", dabs, drel);
//            converged = converged && ((dabs < epsabs) || (drel < epsrel));
//            dabs = max(abs(xi4-xi4_last));
//            drel = max(abs(xi4-xi4_last)/(abs(xi4)+epsabs));
//            printf("  xi_4: dabs = %g, drel = %g\n", dabs, drel);
//            converged = converged && ((dabs < epsabs) || (drel < epsrel));

            xi0_last = xi0;
            xi2_last = xi2;
            xi4_last = xi4;
        }
        else {
            MPI_Reduce(&ndata, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&nrandom, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.dd[0], NULL, nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.dr[0], NULL, nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);
            MPI_Reduce(&binner.rr[0], NULL, nsbins*nmubins, MPI_DOUBLE, MPI_SUM, 0, comm);
        }

        MPI_Bcast(&nconv, 1, MPI_INT, 0, comm);
    }

    if(me == 0) {
        if(nconv < 5)
            printf("Failed to converge after %d samplings.\n", maxsamp);
        printf("Writing estimates to '%s'\n", outfile);
    }

#if 0
    /* Write estimates to file */
    if(me == 0) {
        fprintf(fout, "# Real-space xi(r) estimate\n");
        fprintf(fout, "# Config options: N = %d, M = %d, L = %g, seed = %lu\n", N, M, L, seed);
        fprintf(fout, "# (r values give midpoints of bins)\n");
        fprintf(fout, "# ndd = %g, ndr = %g, nrr = %g\n", estimator.ndd, estimator.ndr, estimator.nrr);
        for(int i = 0; i < estimator.nbins; i++)
            fprintf(fout, "%6.2f %e %e %e\n", estimator.bin_center(i), estimator.simple_xi(i), estimator.dd[i], estimator.rr[i]);
    }
#endif

    if(me == 0) {
        fprintf(fout, "# Plane-parallel xi_ell(s) estimates\n");
        fprintf(fout, "# Config options: N = %d, M = %d, L = %g, seed = %lu\n", N, M, L, seed);
        fprintf(fout, "# (s values give midpoints of bins)\n");
        fprintf(fout, "# ndd = %g, ndr = %g, nrr = %g\n", estimator.ndd, estimator.ndr, estimator.nrr);
        fprintf(fout, "# Columns are: s, xi0, xi2, xi4\n");
        for(int i = 0; i < nsbins; i++)
            fprintf(fout, "%6.2f %e %e %e\n", smin + (i+0.5)*(smax-smin)/nsbins, xi0[i], xi2[i], xi4[i]);
    }

    /* Clean up */
    rho.cleanup();
    vx.cleanup();
    vy.cleanup();
    vz.cleanup();
    if(me == 0)
        fclose(fout);
    MPI_Finalize();

    return 0;
}
