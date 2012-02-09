#ifndef ESTIMATOR_H
#define ESTIMATOR_H

#include "array.h"
#include "particle.h"
#include "vec3.h"

/* Infinity, or at least much much larger than any realistic length scale. */
static const double inf = 1e100;

struct Binner {
    /* Number of pairs in each bin, stored as a double.  Integers less than
     * 2^53 ~= 10^16 have exact representations as doubles, so this only
     * becomes a problem if there are more than 2^53 particles in a single bin. */
    array dd, dr, rr;

    /* Periodic lengths of computational box */
    double Lx, Ly, Lz;

    Binner(double Lx, double Ly, double Lz);
    virtual ~Binner();

    /* (pure virtual function) */
    virtual void add_pair(const Particle& p1, const Particle& p2, array& paircount) = 0;

    /* Correlate pairs of particles from different chaining meshes. */
    virtual void correlate(ChainingMesh& mesh1, ChainingMesh& mesh2, array& paircount, double smax = inf);

    virtual void correlate_dd(ChainingMesh& mesh1, ChainingMesh& mesh2, double smax = inf);
    virtual void correlate_dr(ChainingMesh& mesh1, ChainingMesh& mesh2, double smax = inf);
    virtual void correlate_rr(ChainingMesh& mesh1, ChainingMesh& mesh2, double smax = inf);

    /* Correlate pairs of particles within the same chaining mesh. */
    virtual void self_correlate(ChainingMesh& mesh, array& paircount, double smax = inf);

    virtual void self_correlate_dd(ChainingMesh& mesh, double smax = inf);
    virtual void self_correlate_rr(ChainingMesh& mesh, double smax = inf);
};

/* Bins particle pairs (\vec{s}_1,\vec{s}_2) into bins of s and \mu, assuming
 * the plane-parallel approximation.  The computational box is treated as
 * periodic in the x and y directions, but not in z. */
struct PlaneParallelBinner : public Binner {
    double L;                   // size of computational box

    /* Binning:
     * - s is the redshift-space separation of the pair:
     *     s = |\vec{s}_2 - \vec{s}_1|
     * - mu is the cosine of the angle made by pair separation vector and the
     *   z-axis:
     *     \mu = (s_{2z} - s_{1z})/s */
    int nsbins, nmubins;
    double smin, smax;
    double mumin, mumax;        // always set to 0 and 1

    PlaneParallelBinner(double L, int nsbins, int nmubins, double smin = 0., double smax = 200.);

    virtual void add_pair(const Particle& p1, const Particle& p2, array& paircount);
};


/* Estimates the plane-parallel redshift-space correlation function
 * \xi_s(r,\mu) by counting the number of particle pairs in different
 * configurations. */
struct PlaneParallelEstimator {
    enum {
        PeeblesHauser   = 0,    // DD/RR - 1
        DavisPeebles    = 1,    // DD/DR - 1
        Hamilton        = 2,    // ?
        LandySzalay     = 3     // (DD - 2DR + RR)/RR
    };

    int method;
    int nsbins, nmubins;
    double ndd, ndr, nrr;
    array dd, dr, rr;

    PlaneParallelEstimator(int nsbins, int nmubins, int method = LandySzalay);

    /* Peebles-Hauser DD/RR-1 estimator */
    double ph_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : DD/RR - 1.;
    }

    /* Davis-Peebles DD/DR-1 estimator */
    double dp_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double DR = dr(i1,i2)/ndr;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : (DD - 2*DR + RR)/RR;
    }

    /* Hamilton DD*RR/DR^2-1 estimator */
    double ham_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double DR = dr(i1,i2)/ndr;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : (DD - 2*DR + RR)/RR;
    }

    /* Landy-Szalay (D-R)^2/RR estimator */
    double ls_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double DR = dr(i1,i2)/ndr;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : (DD - 2*DR + RR)/RR;
    }

    /* Compute the Legendre moment \xi_\ell(s) of the correlation function
     * \xi(s,\mu).  The array will have length nsbins. */
    array legendre(int ell);
};


#if 0
/* Estimates the real-space 2-point correlation function \xi(r) by counting the
 * number of particle pairs at different separations. */
struct RealEstimator {
    double L;           // dimensions of computational box

    /* \xi is measured in real-space as a function of r */
    int nbins;          // number of r bins
    double rmin, rmax;  // min and max r

    /* Pair counts, broken into bins.  These are stored as doubles, which
     * ensures an exact representation of integers up to 2^53 or about 10^18.
     * We should never reach a situation where the count within a single bin
     * exceeds this number. */
    array dd, dr, rr;
    double ndd, ndr, nrr;

    void initialize(double L_, int nbins_, double rmin_, double rmax_) {
        L = L_;
        nbins = nbins_;
        rmin = rmin_;
        rmax = rmax_;
        assert(nbins > 0 && rmin < rmax);
        ndd = ndr = nrr = 0;
        dd = array::zeros(nbins);
        dr = array::zeros(nbins);
        rr = array::zeros(nbins);
    }

    void add(const Particle& p1, const Particle& p2, array& count) {
        double r = distance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, L, L, L);

        int i = int(nbins * (r - rmin)/(rmax - rmin));
        if(i >= 0 && i < nbins)
            count[i] += 1;
    }

    void add_dd(const Particle& p1, const Particle& p2) {
        add(p1, p2, dd);
    }

    void add_dr(const Particle& p1, const Particle& p2) {
        add(p1, p2, dr);
    }

    void add_rr(const Particle& p1, const Particle& p2) {
        add(p1, p2, rr);
    }


    /* Correlate pairs of particles within different ChainingMeshes.  Particle
     * type specifier is 0 for data, 1 for random. */
    void correlate(ChainingMesh& mesh1, int t1, ChainingMesh& mesh2, int t2, double vmax = 0.) {
        /* Use self_correlate() for correlating a ChainingMesh of particles
         * with itself. */
        assert(&mesh1 != &mesh2);

        /* Choose which count these particles contribute to (DD, DR, or RR) */
        array* count;
        if(!t1 && !t2)
            count = &dd;
        else if(t1 && t2)
            count = &rr;
        else
            count = &dr;

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
                 * r that we actually care about */
                if(cell_separation(*c1, *c2, L, L, L) > rmax)
                    continue;

                /* Correlate particles within these cells */
                for(Particle* p1 = p1begin; p1 < p1end; p1++)
                    for(Particle* p2 = p2begin; p2 < p2end; p2++)
                        add(*p1, *p2, *count);
            }
        }
    }

    /* Count pairs of particles within the same ChainingMesh.  Particle type
     * specifier is 0 for data, 1 for random. */
    void self_correlate(ChainingMesh& mesh, int t, double vmax = 0.) {
        array* count;
        if(!t)
            count = &dd;
        else
            count = &rr;

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
                if(cell_separation(*c1, *c2, L, L, L) > rmax)
                    continue;

                /* Correlate particles within these cells */
                for(Particle* p1 = p1begin; p1 < p1end; p1++)
                    for(Particle* p2 = p2begin; p2 < p2end; p2++)
                        add(*p1, *p2, *count);
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
                    add(*p1, *p2, *count);
        }
    }

    /* Simple DD/RR estimator */
    double simple_xi(int i) {
        assert(i >= 0 && i < nbins);
        if(rr[i] == 0)
            return 0.;
        else
            return (dd[i]/ndd) / (rr[i]/nrr) - 1;
    }

    double bin_center(int i) {
        assert(i >= 0 && i < nbins);
        return rmin + (i + 0.5)*(rmax - rmin)/nbins;
    }
};
#endif // 0


#endif // ESTIMATOR_H
