#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "array.h"
#include "particle.h"
#include "vec3.h"


/* Estimates the redshift-space correlation function \xi_s(r,\mu) within the
 * plane-parallel approximation, by counting the number of particle pairs in
 * different configurations. */
struct StandardEstimator {
    double L;                   // size of computational box

    /* \xi_s is treated as a function of s and mu:
     * - s is the redshift-space separation of the pair
     * - mu is the cosine of the angle made by pair separation vector and the
     *   z-axis */
    int nsbins, nmubins;
    double smin, smax;
    double mumin, mumax;

    /* Number of pairs in each (s,mu)-bin, stored as a double.  Integers less
     * than 2^53 ~= 10^16 have exact representations as doubles, so this only
     * becomes a problem if there are more than 2^53 particles in a single bin. */
    array dd, dr, rr;

    /* The total number of pairs of each type.  These values must be set by
     * hand, since geometrical culling (e.g. a chaining mesh) will cause a
     * large number of pairs to be skipped over without being counted
     * explicitly. */
    double ndd, ndr, nrr;

    StandardEstimator(double L,
                      int nsbins, double smin, double smax,
                      int nmubins, double mumin = 0., double mumax = 1.);

    void add_pair(const Particle& p1, const Particle& p2, array& paircount);

    /* Correlate pairs of particles from different chaining meshes. */
    void correlate(ChainingMesh& mesh1, ChainingMesh& mesh2, array& paircount, double vmax = 0.);
    void correlate_dd(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax = 0.);
    void correlate_dr(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax = 0.);
    void correlate_rr(ChainingMesh& mesh1, ChainingMesh& mesh2, double vmax = 0.);

    /* Correlate pairs of particles within the same chaining mesh. */
    void self_correlate(ChainingMesh& mesh, array& paircount, double vmax = 0.);
    void self_correlate_dd(ChainingMesh& mesh, double vmax = 0.);
    void self_correlate_rr(ChainingMesh& mesh, double vmax = 0.);

    /* Peebles-Hauser DD/RR estimator */
    double ph_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : DD/RR - 1.;
    }

    /* Landy-Szalay (D-R)^2/RR estimator */
    double ls_xi(int i1, int i2) {
        double DD = dd(i1,i2)/ndd;
        double DR = dr(i1,i2)/ndr;
        double RR = rr(i1,i2)/nrr;
        return (RR == 0) ? 0. : (DD - 2*DR + RR)/RR;
    }

    /* Compute the Legendre moment \xi_\ell(s) of the correlation function
     * \xi(s,\mu).  Note that this calculation only really makes sense if
     * mumin = 0 and mumax = 1. */
    void legendre(int ell, array& s, array& xi);

    void bin_center(int i1, int i2, double& s, double& mu) {
        assert(i1 >= 0 && i1 < nsbins && i2 >= 0 && i2 < nmubins);
        s = smin + (i1 + 0.5)*(smax - smin)/nsbins;
        mu = mumin + (i2 + 0.5)*(mumax - mumin)/nmubins;
    }
};

#if 0
/* Estimates the 2-point correlation function by counting the number of
 * particle pairs in different configurations. */
struct ParallelRedshiftXiEstimator {
    double L;                   // size of computational box

    /* \xi is treated as a function of s and mu:
     * - s is the redshift-space pair separation of the pair
     * - mu is the cosine of the angle made by pair separation vector and the
     *   z-axis */
    int nsbins, nmubins;
    double smin, smax;
    double mumin, mumax;

    /* The number of pairs in each configuration is stored as a double (ints
     * less than 2^53 ~ 10^16 have exact representations as doubles). */
    array dd, dr, rr;
    double ndd, ndr, nrr;

    void initialize(double L_,
                    int nsbins_, double smin_, double smax_,
                    int nmubins_, double mumin_, double mumax_)
    {
        L = L_;
        nsbins = nsbins_;
        smin = smin_;
        smax = smax_;
        nmubins = nmubins_;
        mumin = mumin_;
        mumax = mumax_;
        assert(nsbins > 0 && nmubins > 0 && smin < smax && mumin < mumax);
        ndd = ndr = nrr = 0;
        dd = array::zeros(nsbins, nmubins);
        dr = array::zeros(nsbins, nmubins);
        rr = array::zeros(nsbins, nmubins);
    }

    void add(const Particle& p1, const Particle& p2, array& count) {
        vec3d s1(p1.x, p1.y, p1.z + p1.vz);
        vec3d s2(p2.x, p2.y, p2.z + p2.vz);
        vec3d S = s2 - s1;

        double s = distance(p1.x, p1.y, p1.z + p1.vz, p2.x, p2.y, p2.z + p2.vz, L, L);
        double mu = fabs(p2.z + p2.vz - p1.z - p1.vz)/s;

        int i1 = int(nsbins * (s - smin)/(smax - smin));
        int i2 = int(nmubins * (mu - mumin)/(mumax - mumin));
        if(i1 >= 0 && i1 < nsbins && i2 >= 0 && i2 < nmubins)
            count(i1, i2) += 1;
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
                if(cell_separation(*c1, *c2, L, L) > smax + 2*vmax)
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
                if(cell_separation(*c1, *c2, L, L) > smax + 2*vmax)
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
    double simple_xi(int i1, int i2) {
        if(rr(i1,i2) == 0)
            return 0.;

        double DD = dd(i1,i2)/ndd;
        double RR = rr(i1,i2)/nrr;
        return DD/RR - 1.;
    }

    /* Landy-Szalay estimator */
    double ls_xi(int i1, int i2) {
        if(rr(i1,i2) == 0)
            return 0.;

        double DD = dd(i1,i2)/ndd;
        double DR = dr(i1,i2)/ndr;
        double RR = rr(i1,i2)/nrr;
        return (DD - 2*DR + RR)/RR;
    }

    /* Compute the Legendre moment \xi_\ell(s) of the correlation function
     * \xi(s,mu).  Note that this calculation only really makes sense if
     * mumin = 0 and mumax = 1. */
    void legendre(int ell, array& s, array& xi) {
        s.resize(nsbins);
        xi.resize(nsbins);
        double mu, plmu;
        for(int i1 = 0; i1 < nsbins; i1++) {
            s[i1] = smin + (i1 + 0.5)*(smax - smin)/nsbins;
            xi[i1] = 0.;
            for(int i2 = 0; i2 < nmubins; i2++) {
                mu = mumin + (i2 + 0.5)*(mumax - mumin)/nmubins;
                plmu = LegendrePolynomial(ell, mu);
                xi[i1] += simple_xi(i1, i2)*plmu;
            }
            xi[i1] *= (2*ell+1) * (mumax - mumin)/nmubins;
        }
    }

    void bin_center(int i1, int i2, double& s, double& mu) {
        assert(i1 >= 0 && i1 < nsbins && i2 >= 0 && i2 < nmubins);
        s = smin + (i1 + 0.5)*(smax - smin)/nsbins;
        mu = mumin + (i2 + 0.5)*(mumax - mumin)/nmubins;
    }
};

/* Estimates the real-space 2-point correlation function \xi(r) by counting the
 * number of particle pairs at different separations. */
struct RealXiEstimator {
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


#endif // HISTOGRAM_H
