#include <cassert>
#include <cmath>

#include "estimator.h"
#include "vec3.h"


/* Return the distance between the points (x1,y1,z1) and (x2,y2,z2), assuming
 * the points lie within a periodic box with side lengths Lx, Ly, Lz. */
static double distance(double x1, double y1, double z1,
                       double x2, double y2, double z2,
                       double Lx = inf, double Ly = inf, double Lz = inf)
{
    double rx = fabs(x2 - x1);
    double ry = fabs(y2 - y1);
    double rz = fabs(z2 - z1);
    rx = fmin(rx, Lx - rx);
    ry = fmin(ry, Ly - ry);
    rz = fmin(rz, Lz - rz);
    return sqrt(rx*rx + ry*ry + rz*rz);
}

/* Find the smallest distance between intervals (a1,b1) and (a2,b2), which lie
 * within the periodic space [0,L]. */
static double interval_separation(double a1, double b1, double a2, double b2, double L = inf) {
    assert(a1 <= b1 && a2 <= b2);

    if(a1 <= a2)
        return fmax(0., fmin(a2 - b1, L - b2 + a1));
    else
        return fmax(0., fmin(a1 - b2, L - b1 + a2));
}

/* Find the smallest distance between two cells, within the periodic box with
 * side lengths Lx, Ly, Lz. */
static double cell_separation(const Cell& c1, const Cell& c2, double Lx = inf, double Ly = inf, double Lz = inf) {
    double dx = interval_separation(c1.xmin, c1.xmax, c2.xmin, c2.xmax, Lx);
    double dy = interval_separation(c1.ymin, c1.ymax, c2.ymin, c2.ymax, Ly);
    double dz = interval_separation(c1.zmin, c1.zmax, c2.zmin, c2.zmax, Lz);
    return sqrt(dx*dx + dy*dy + dz*dz);
}


/***** Binner *****/

Binner::Binner(double Lx_, double Ly_, double Lz_) {
    Lx = Lx_;
    Ly = Ly_;
    Lz = Lz_;
}

Binner::~Binner() {
}

void Binner::correlate(ChainingMesh& mesh1, ChainingMesh& mesh2, array& paircount, double smax) {
    /* Iterate over cells in mesh 1... */
    for(vector<Cell>::iterator c1 = mesh1.cells.begin(); c1 != mesh1.cells.end(); c1++) {
        Particle* p1begin = &mesh1.particles[c1->start];
        Particle* p1end = p1begin + c1->count;
        if(p1begin == p1end)
            continue;

        /* ... and cells in mesh 2 */
        for(vector<Cell>::iterator c2 = mesh2.cells.begin(); c2 != mesh2.cells.end(); c2++) {
            Particle* p2begin = &mesh2.particles[c2->start];
            Particle* p2end = p2begin + c2->count;
            if(p2begin == p2end)
                continue;

            /* Check distance between cells, and compare it to the maximum
             * separation that we actually care about */
            if(cell_separation(*c1, *c2, Lx, Ly, Lz) > smax)
                continue;

            /* Correlate particles within these cells */
            for(Particle* p1 = p1begin; p1 < p1end; p1++)
                for(Particle* p2 = p2begin; p2 < p2end; p2++)
                    add_pair(*p1, *p2, paircount);
        }
    }
}

void Binner::self_correlate(ChainingMesh& mesh, array& paircount, double smax) {
    /* Iterate over pairs of distinct cells within the ChainingMesh */
    for(vector<Cell>::iterator c1 = mesh.cells.begin(); c1 != mesh.cells.end(); c1++) {
        Particle* p1begin = &mesh.particles[c1->start];
        Particle* p1end = p1begin + c1->count;
        if(p1begin == p1end)
            continue;

        for(vector<Cell>::iterator c2 = c1+1; c2 != mesh.cells.end(); c2++) {
            Particle* p2begin = &mesh.particles[c2->start];
            Particle* p2end = p2begin + c2->count;
            if(p2begin == p2end)
                continue;

            /* Check distance between cells, and compare it to the maximum
             * r that we actually care about */
            if(cell_separation(*c1, *c2, Lx, Ly, Lz) > smax)
                continue;

            /* Correlate particles within these cells */
            for(Particle* p1 = p1begin; p1 < p1end; p1++)
                for(Particle* p2 = p2begin; p2 < p2end; p2++)
                    add_pair(*p1, *p2, paircount);
        }
    }

    /* Now count pairs of distinct particles from within the same cell */
    for(vector<Cell>::iterator c = mesh.cells.begin(); c != mesh.cells.end(); c++) {
        Particle* pbegin = &mesh.particles[c->start];
        Particle* pend = pbegin + c->count;
        if(pbegin == pend)
            continue;

        for(Particle* p1 = pbegin; p1 != pend; p1++)
            for(Particle* p2 = p1+1; p2 != pend; p2++)
                add_pair(*p1, *p2, paircount);
    }
}

void Binner::correlate_dd(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax) {
    correlate(mesh1, mesh2, dd, vmax);
}

void Binner::correlate_dr(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax) {
    correlate(mesh1, mesh2, dr, vmax);
}

void Binner::correlate_rr(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax) {
    correlate(mesh1, mesh2, rr, vmax);
}

void Binner::self_correlate_dd(ChainingMesh& mesh, double vmax) {
    self_correlate(mesh, dd, vmax);
}

void Binner::self_correlate_rr(ChainingMesh& mesh, double vmax) {
    self_correlate(mesh, rr, vmax);
}


/***** PlaneParallelBinner *****/

PlaneParallelBinner::PlaneParallelBinner
    (double L, int nsbins_, int nmubins_, double smin_, double smax_)
    : Binner(L, L, inf)
{
    nsbins = nsbins_;
    nmubins = nmubins_;
    smin = smin_;
    smax = smax_;
    mumin = 0.;
    mumax = 1;
    assert(nsbins > 0 && nmubins > 0 && smin < smax);
    dd = array::zeros(nsbins, nmubins);
    dr = array::zeros(nsbins, nmubins);
    rr = array::zeros(nsbins, nmubins);
}

void PlaneParallelBinner::add_pair(const Particle& p1, const Particle& p2, array& paircount) {
    double s = distance(p1.x, p1.y, p1.z + p1.vz,
                        p2.x, p2.y, p2.z + p2.vz,
                        Lx, Ly, Lz);
    double mu = fabs(p2.z + p2.vz - p1.z - p1.vz)/s;    // use symmetry between +mu and -mu

    int i1 = int(nsbins * (s - smin)/(smax - smin));
    int i2 = int(nmubins * (mu - mumin)/(mumax - mumin));
    if(i1 >= 0 && i1 < nsbins && i2 >= 0 && i2 < nmubins)
        paircount(i1, i2) += 1;
}


/***** PlaneParallelEstimator *****/

PlaneParallelEstimator::PlaneParallelEstimator(int nsbins_, int nmubins_, int method_) {
    nsbins = nsbins_;
    nmubins = nmubins_;
    method = method_;
    ndd = ndr = nrr = 0;
    dd = array::zeros(nsbins, nmubins);
    dr = array::zeros(nsbins, nmubins);
    rr = array::zeros(nsbins, nmubins);
}

array PlaneParallelEstimator::legendre(int ell) {
    /* Calculate Legendre polynomial values */
    double mumin = 0., mumax = 1.;
    array plmu(ell+1, nmubins);
    if(ell == 0)
        plmu = array::ones(1, nmubins);
    else {
        double x;
        for(int j = 0; j < nmubins; j++) {
            x = mumin + (j + 0.5)*(mumax - mumin)/nmubins;
            plmu(0,j) = 1.;
            plmu(1,j) = x;
            for(int n = 1; n < ell; n++)
                plmu(n+1,j) = ((2*n+1)*x*plmu(n,j) - n*plmu(n-1,j))/(n+1);
        }
    }

    array xi = array::zeros(nsbins);
    for(int i = 0; i < nsbins; i++) {
        switch(method) {
            case PeeblesHauser:
                for(int j = 0; j < nmubins; j++)
                    xi[i] += ph_xi(i, j) * plmu(ell,j);
                break;
            case LandySzalay:
                for(int j = 0; j < nmubins; j++)
                    xi[i] += ls_xi(i, j) * plmu(ell,j);
                break;
            default:
                break;
        }

        xi[i] *= (2*ell+1) * (mumax - mumin)/nmubins;
    }
    return xi;
}
