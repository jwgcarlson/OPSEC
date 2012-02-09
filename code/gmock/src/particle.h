#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
using std::vector;

/* Particle position in real or redshift space */
struct Particle {
    float x, y, z;
    float vx, vy, vz;
};

/* A cell within the chaining mesh.  Note that this Cell is implicitly tied to
 * a ChainingMesh, even though there is no explicit reference. */
struct Cell {
    /* Physical extent of cell in real space */
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;

    /* Starting index of particles for this cell */
    int start;

    /* Number of particles in this cell */
    int count;
};

/* A set of particles within a slab of the computational box, arranged into a
 * chaining mesh. */
struct ChainingMesh {
    vector<Cell> cells;         // list of all mesh cells within slab, indexed by j = (jx*my + jy)*mz + jz
    vector<Particle> particles; // list of all particles within slab, ordered my cell
};

#endif // PARTICLE_H
