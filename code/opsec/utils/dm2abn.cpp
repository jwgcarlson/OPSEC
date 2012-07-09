/* dm2abn
 *
 * Takes a dm_X.bin file of subsampled simulation particles and constructs a
 * galaxy sample of the form required by KLPK.  A dm_X.bin file is a binary
 * format consisting of a header
 *   int   endian;         // Endianness flag, should equal 1 (if endianness matches)
 *   int   hsize;          // Size of remaining header, should equal 4*sizeof(int) + 2*sizeof(float)
 *   int   npart;          // Total number of particles.
 *   int   nsph;           // Number of gas particles.  
 *   int   nstar;          // Number of star particles. 
 *   float aa;             // Scale factor.
 *   float softlen;        // Gravitational softening.
 * followed by 3*npart floats giving the positions of the particles, 3*npart
 * floats giving their velocities, and npart floats/ints specifying some other
 * scalar property.  Positions are specified in scaled box coordinates, i.e.
 * with 0 <= x,y,z < 1, so that the size of the simulation box must be known in
 * order to construct physical positions.  Velocities are given in units of
 * peculiar velocity over aH.  (XXX This may be totally false, I don't know how
 * the velocities are stored, but whatever it is it should be documented here.) */

#include <cmath>
#include <cstdio>
#include <cstring>
#include <unistd.h>

#include <string>
#include <vector>
using std::string;
using std::vector;

#include "abn.h"
#include "cfg.h"
#include "vec3.h"

/* Struct representing a galaxy either in Cartesian or spherical coordinates */
struct Galaxy {
    union { float x; float r; };
    union { float y; float mu; };
    union { float z; float phi; };
};

static inline bool contains(double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax, const Galaxy& g) {
    return (XMin <= g.x && g.x <= XMax)
        && (YMin <= g.y && g.y <= YMax)
        && (ZMin <= g.z && g.z <= ZMax);
}

/* Construct galaxy list in spherical coordinates */
void read_galaxies_spherical(int npart, FILE* fpos, FILE* fvel, float L,
                             double RMin, double RMax, double MuMin, double MuMax, double PhiMin, double PhiMax,
                             bool rsd, const vec3f& origin,
                             vector<Galaxy>& galaxies)
{
    int i = 0;
    int nread;
    vec3f pos, vel, n;
    Galaxy g;
    while(i < npart && !feof(fpos) && !feof(fvel)) {
        /* Read position and velocity of particle */
        nread = fread(&pos, sizeof(float), 3, fpos);
        nread = fread(&vel, sizeof(float), 3, fvel);

        /* Rescale and translate */
        pos = L*pos - origin;

        /* Apply redshift distortion */
        if(rsd) {
            n = unit(pos);
            pos += n * dot(n,vel) * L;
        }

        /* Compute spherical coordinates */
        g.r = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
        g.mu = pos.z/g.r;
        g.phi = fmod(M_PI + atan2(pos.y, pos.x), 2*M_PI);

        if(contains(RMin, RMax, MuMin, MuMax, PhiMin, PhiMax, g))
            galaxies.push_back(g);
        i++;
    }

    if(i != npart)
        fprintf(stderr, "read_galaxies_spherical: error encountered while reading particles (%d/%d read)\n", i+1, npart);
}

/* Construct galaxy list in Cartesian coordinates */
/* XXX I'm not sure what I want this method to do.  Should it assume the
 * distant-observer approximation, and apply redshifts only along one axis?  Or
 * should it do exactly the same as the spherical coordinate case, but with
 * final galaxy positions in Cartesian coordinates?  I feel like the only
 * unambiguous use case is in the absence of redshift-space distortions, so
 * should I make this the only functionality?  For now, I'm taking this route,
 * but it should definitely be revisited later. */
void read_galaxies_cartesian(int npart, FILE* fpos, FILE* fvel, float L,
                             double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax,
                             bool rsd, const vec3f& origin,
                             vector<Galaxy>& galaxies)
{
    if(rsd)
        fprintf(stderr, "read_galaxies_cartesian: warning: RSD not implemented for Cartesian coordinates\n");

    int i = 0;
    int nread;
    vec3f pos, n;
    Galaxy g;
    while(i < npart && !feof(fpos) && !feof(fvel)) {
        /* Read position of particle */
        nread = fread(&pos, sizeof(float), 3, fpos);

        /* Rescale and translate */
        pos = L*pos;

        /* Compute spherical coordinates */
        g.x = pos.x;
        g.y = pos.y;
        g.z = pos.z;

        if(contains(XMin, XMax, YMin, YMax, ZMin, ZMax, g))
            galaxies.push_back(g);

        i++;
    }
}

int main(int argc, char* argv[]) {
    Config cfg = cfg_new();

    /* Parse command line switches */
    int opt;
    const char* optstring = "hc:";
    while((opt = getopt(argc, argv, optstring)) != -1) {
        switch(opt) {
        case 'h':
            printf("Usage: %s [-c <config>] OPTIONS\n", argv[0]);
            return 0;
        case 'c':
            cfg_read_file(cfg, optarg);
            break;
        default:
            printf("Usage: %s [-c config] OPTIONS\n", argv[0]);
            return 1;
        }
    }

    /* Parse additional command line options */
    for(int i = optind; i < argc; i++)
        cfg_read_line(cfg, argv[i]);

    /* Make sure output file is specified and that we can write to it */
    if(!cfg_has_key(cfg, "galfile")) {
        fprintf(stderr, "dm2abn: galfile not specified in configuration\n");
        return 1;
    }
    const char* galfile = cfg_get(cfg, "galfile");
    FILE* fout = fopen(galfile, "w");
    if(!fout) {
        fprintf(stderr, "dm2abn: could not open '%s' for writing.\n", galfile);
        return 1;
    }

    /* Check coordinate system */
    string coordsys(cfg_get(cfg, "coordsys"));
    if(coordsys != "spherical" && coordsys != "cartesian") {
        fprintf(stderr, "dm2abn: missing or invalid coordsys option '%s'\n", coordsys.c_str());
        return 1;
    }

    /* Read remaining config options */
    if(char* missing = cfg_missing_keys(cfg, "dmfile,Lbox,origin")) {
        opsec_error("dm2abn: missing config options %s\n", missing);
        return 1;
    }
    const char* dmfile = cfg_get(cfg, "dmfile");
    double Lbox = cfg_get_double(cfg, "Lbox");
    vec3f origin;
    if(3 != sscanf(cfg_get(cfg, "origin"), "%f %f %f", &origin.x, &origin.y, &origin.z)) {
        opsec_error("dm2abn: invalid value for config option 'origin': %s\n", cfg_get(cfg, "origin"));
        return 1;
    }

    /* Open dm file */
    FILE* fdm;
    fdm = fopen(dmfile, "r");
    if(fdm == NULL) {
        opsec_error("dm2abn: could not open file '%s'\n", dmfile);
        return 0;
    }

    /* Read endian flag */
    int endian;
    int nread = fread(&endian, sizeof(int), 1, fdm);
    if(endian != 1) {
        opsec_error("dm2abn: endian flag is not 1\n");
        fclose(fdm);
        return 0;
    }

    struct FileHeader {
        int   npart;          /* Total number of particles. */
        int   nsph;           /* Number of gas particles.   */
        int   nstar;          /* Number of star particles.  */
        float aa;             /* Scale factor. */
        float softlen;        /* Gravitational softening    */
    };

    /* Read file header */
    int hsize;
    nread = fread(&hsize, sizeof(int), 1, fdm);
    if(nread != 1 || hsize != sizeof(struct FileHeader)) {
        opsec_error("dm2abn: header size is %d, expecting %zd\n", hsize, sizeof(struct FileHeader));
        fclose(fdm);
        return 0;
    }
    FileHeader header;
    nread = fread(&header, hsize, 1, fdm);

    /* Helpful information */
    printf("%d dark matter particles in a %g^3 box: nbar = %g\n", header.npart, Lbox, header.npart/(Lbox*Lbox*Lbox));

    /* Open new handle for reading particle velocities */
    long start = ftell(fdm);
    FILE* fvel = fopen(dmfile, "r");
    fseek(fvel, start + 3*header.npart*sizeof(float), SEEK_SET);

    /* No redshift distortions yet */
    bool rsd = false;

    vector<Galaxy> galaxies;
    if(coordsys == "spherical") {
        if(char* missing = cfg_missing_keys(cfg, "RMin,RMax,MuMin,MuMax,PhiMin,PhiMax")) {
            fprintf(stderr, "dm2abn: missing config options %s\n", missing);
            return 1;
        }
        double RMin = cfg_get_double(cfg, "RMin");
        double RMax = cfg_get_double(cfg, "RMax");
        double MuMin = cfg_get_double(cfg, "MuMin");
        double MuMax = cfg_get_double(cfg, "MuMax");
        double PhiMin = cfg_get_double(cfg, "PhiMin");
        double PhiMax = cfg_get_double(cfg, "PhiMax");
        read_galaxies_spherical(header.npart, fdm, fvel, Lbox, RMin, RMax, MuMin, MuMax, PhiMin, PhiMax, rsd, origin, galaxies);
    }
    else if(coordsys == "cartesian") {
        if(!cfg_has_keys(cfg, "XMin,XMax,YMin,YMax,ZMin,ZMax", ",")) {
            fprintf(stderr, "dm2abn: must specify config options {X,Y,Z}{Min,Max}\n");
            return 1;
        }
        double XMin = cfg_get_double(cfg, "XMin");
        double XMax = cfg_get_double(cfg, "XMax");
        double YMin = cfg_get_double(cfg, "YMin");
        double YMax = cfg_get_double(cfg, "YMax");
        double ZMin = cfg_get_double(cfg, "ZMin");
        double ZMax = cfg_get_double(cfg, "ZMax");
        read_galaxies_cartesian(header.npart, fdm, fvel, Lbox, XMin, XMax, YMin, YMax, ZMin, ZMax, rsd, origin, galaxies);
    }

    fclose(fdm);
    fclose(fvel);

    /* Write galaxies to file */
    Config opts = cfg_new();
    cfg_set(opts, "coordsys", coordsys.c_str());
    fprintf(fout, "# Galaxies converted to .abn format from '%s'\n", dmfile);
    abn_write(fout, &galaxies[0], galaxies.size(), "3f", opts);
    fclose(fout);

    return 0;
}
