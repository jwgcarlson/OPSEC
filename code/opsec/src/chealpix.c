/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *  Refactoring and cleanup by Jordan Carlson <jwgcarlson@gmail.com>.
 *  Based on version 2.14a of the official HEALPix software distribution.
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */

/* Standard Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local Includes */
#include "chealpix.h"
#include "fitsio.h"

const double pi = M_PI;
const double twopi = 2*M_PI;
const double piover2 = 0.5*M_PI;

#define MAX_NSIDE 8192

/* Macro to print file and line number in addition to error message */
//#define eprintf_helper(...) fprintf(stderr, __VA_ARGS__)
#define eprintf(fmt, ...) fprintf(stderr, "%s (%d): " fmt, __FILE__, __LINE__, ##__VA_ARGS__)

/* Sets the array giving the number of the pixel lying in (x,y).
 * - x and y are in {1,128}
 * - the pixel number is in {0,128**2-1} */
void mk_xy2pix(int *x2pix, int *y2pix) {
  int id, ip, i, j, k;
  for(i = 0; i < 128; i++) {
    j = i;
    k = 0;
    ip = 1;
    while(j > 0) {
      id = j % 2;
      j /= 2;
      k += ip*id;
      ip *= 4;
    }
    x2pix[i] = k;
    y2pix[i] = 2*k;
  }
}

/* Constructs arrays giving x and y in the face from pixel number for the
 * nested (quad-cube like) ordering of pixels.  The bits corresponding to x and
 * y are interleaved in the pixel number.  One breaks up the pixel number by
 * even and odd bits. */
void mk_pix2xy(int *pix2x, int *pix2y) {
  int id, ip, ix, iy, j, k;
  for(k = 0; k < 1024; k++) {
    j = k;
    ix = iy = 0;
    ip = 1;         /* bit position (in x and y) */
    while(j != 0) { /* go through all the bits */
      id = j % 2;   /* bit value (in k), goes in ix */
      j /= 2;
      ix += id*ip;
      id = j % 2;
      j /= 2;
      iy += id*ip;
      ip *= 2;
    }
    pix2x[k] = ix;    /* in 0,31 */
    pix2y[k] = iy;    /* in 0,31 */
  }
}

/* Gives the pixel number ipix (NESTED) corresponding to angles theta and phi.
 * The computation is made to the highest resolution available (nside=8192) and
 * then degraded to that required (by integer division).  This doesn't cost
 * more, and it makes sure that the treatement of round-off will be consistent
 * for every resolution. */
void ang2pix_nest(long nside, double theta, double phi, long *ipix) {
  double z, za, z0, tt, tp, tmp;
  int face_num,jp,jm;
  long ifp, ifm;
  int ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if(nside < 1 || nside > MAX_NSIDE) {
    eprintf("nside out of range: %ld\n", nside);
    exit(0);
  }
  if(theta < 0. || theta > pi) {
    eprintf("theta out of range: %f\n", theta);
    exit(0);
  }
  if(!setup_done) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if(phi >= twopi)  phi = phi - twopi;
  else if(phi < 0.) phi = phi + twopi;
  tt = phi / piover2; /* in [0,4) */
  
  if(za <= z0) { /* equatorial region */
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(MAX_NSIDE*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(MAX_NSIDE*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / MAX_NSIDE; /* in {0,4} */
    ifm = jm / MAX_NSIDE;
    
    if(ifp == ifm)
      face_num = (ifp % 4) + 4; /* faces 4 to 7 */
    else if(ifp < ifm)
      face_num = ifp % 4;       /* (half-)faces 0 to 3 */
    else
      face_num = (ifm % 4) + 8; /* (half-)faces 8 to 11 */
    
    ix = jm % MAX_NSIDE;
    iy = MAX_NSIDE - (jp % MAX_NSIDE) - 1;
  }
  else { /* polar region, za > 2/3 */
    ntt = (int)floor(tt);
    if(ntt >= 4) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt(3.*(1. - za)); /* in (0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up) */
    /* line going toward the pole as phi increases */
    jp = (int)floor(MAX_NSIDE * tp * tmp); 

    /* that one goes away of the closest pole */
    jm = (int)floor(MAX_NSIDE * (1. - tp) * tmp);
    jp = (int)(jp < MAX_NSIDE-1 ? jp : MAX_NSIDE-1);
    jm = (int)(jm < MAX_NSIDE-1 ? jm : MAX_NSIDE-1);
    
    /* finds the face and pixel's (x,y) */
    if(z >= 0) {
      face_num = ntt; /* in {0,3} */
      ix = MAX_NSIDE - jm - 1;
      iy = MAX_NSIDE - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix = jp;
      iy = jm;
    }
  }
  
  ix_low = ix % 128;
  ix_hi  = ix/128;
  iy_low = iy % 128;
  iy_hi  = iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128) + (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(MAX_NSIDE/nside,2));      /* in {0, nside**2 - 1} */
  *ipix =(long)(ipf + face_num*pow(nside,2));   /* in {0, 12*nside**2 - 1} */
}

/* Gives the pixel number ipix (RING) corresponding to angles theta and phi. */
void ang2pix_ring(long nside, double theta, double phi, long *ipix) {
  int nl2, nl4, ncap, npix, jp, jm, ipix1;
  double z, za, tt, tp, tmp;
  int ir, ip, kshift;
  
  double z0 = 2.0/3.0;
  
  if(nside < 1 || nside > MAX_NSIDE) {
    eprintf("nside out of range: %ld\n", nside);
    exit(0);
  }
  
  if(theta < 0. || theta > pi) {
    eprintf("theta out of range: %f\n", theta);
    exit(0);
  }
  
  z = cos(theta);
  za = fabs(z);
  if(phi >= twopi)  phi -= twopi;
  else if(phi < 0.) phi += twopi;
  tt = phi / piover2;   /* in [0,4) */
  
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = nl2*(nside-1); /* number of pixels in the north polar cap */
  npix = 12*nside*nside;
  
  if(za <= z0) {
    jp = (int)floor(nside*(0.5 + tt - z*0.75)); /* index of ascending edge line */
    jm = (int)floor(nside*(0.5 + tt + z*0.75)); /* index of descending edge line */
    
    ir = nside + 1 + jp - jm;   /* in {1,2n+1} (ring number counted from z=2/3) */
    kshift = 0;
    if((ir % 2) == 0.) kshift = 1;      /* kshift=1 if ir even, 0 otherwise */
    
    ip = (jp+jm - nside + kshift + 1)/2 + 1;    /* in {1,4n} */
    if(ip > nl4) ip = ip - nl4;
    
    ipix1 = ncap + nl4*(ir-1) + ip;
  }
  else {
    tp = tt - floor(tt);    /* MOD(tt,1.d0) */
    tmp = sqrt(3.*(1. - za));
    
    jp = (int)floor(nside * tp * tmp);          /* increasing edge line index */
    jm = (int)floor(nside * (1. - tp) * tmp);   /* decreasing edge line index */
    
    ir = jp + jm + 1;           /* ring number counted from the closest pole */
    ip = (int)floor(tt*ir) + 1; /* in {1,4*ir} */
    if(ip > 4*ir)
        ip -= 4*ir;
    
    ipix1 = 2*ir*(ir-1) + ip;
    if(z <= 0.)
      ipix1 = npix - 2*ir*(ir+1) + ip;
  }
  *ipix = ipix1 - 1;    /* in {0, npix-1} */
}

/* Gives theta and phi corresponding to pixel ipix (NESTED) for a parameter nside. */
void pix2ang_nest(long nside, long ipix, double *theta, double *phi) {
  int npix, npface, face_num;
  int ipf, ip_low, ip_trunc, ip_med, ip_hi;
  int ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
  double z, fn, fact1, fact2;
  
  static char setup_done = 0;
  static int pix2x[1024], pix2y[1024];
  
  /* Coordinate of the lowest corner of each face */
  const int jrll[12] = { 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };  /* in units of nside */
  const int jpll[12] = { 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };  /* in units of nside/2 */

  if(nside < 1 || nside > MAX_NSIDE) {
    eprintf("nside out of range: %ld\n", nside);
    exit(0);
  }
  npix = 12 * nside*nside;
  if(ipix < 0 || ipix > npix-1) {
    eprintf("ipix out of range: %ld\n", ipix);
    exit(0);
  }

  /* Initiate the array for the pixel number -> (x,y) mapping */
  if(!setup_done) {
    mk_pix2xy(pix2x,pix2y);
    setup_done = 1;
  }

  fn = 1.*nside;
  fact1 = 1./(3.*fn*fn);
  fact2 = 2./(3.*fn);
  nl4 = 4*nside;

  /* Find the face, and the number in the face */
  npface = nside*nside;
  face_num = ipix/npface;       /* face number, in {0,11} */
  ipf = ipix % npface;          /* pixel number in the face, in {0,npface-1} */

  /* Find x,y on the face (starting from the lowest corner) from the pxel number */
  ip_low = ipf % 1024;          /* content of the last 10 bits */
  ip_trunc = ipf/1024;          /* truncation of the last 10 bits */
  ip_med = ip_trunc % 1024;     /* content of the next 10 bits */
  ip_hi = ip_trunc/1024;        /* content of the high weight 10 bits */

  ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
  iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

  /* Transform this in (horizontal, vertical) coordinates */
  jrt = ix + iy;        /* 'vertical' in {0,2*(nside-1)} */
  jpt = ix - iy;        /* 'horizontal' in {-nside+1,nside-1} */

  /* Compute the z coordinate on the sphere */
  jr = jrll[face_num]*nside - jrt - 1;  /* ring number in {1,4*nside-1} */
  if(jr < nside) {              /* north pole region */
     nr = jr;
     z = 1. - nr*nr*fact1;
     kshift = 0;
  }
  else if(jr > 3*nside) {       /* south pole region */
     nr = nl4 - jr;
     z = - 1. + nr*nr*fact1;
     kshift = 0;
  }
  else {
    nr = nside;                 /* equatorial region (the most frequent) */
    z = (2*nside-jr)*fact2;
    kshift = (jr - nside) % 2;
  }
  *theta = acos(z);
  
  /* Compute the phi coordinate on the sphere, in [0,2*pi] */
  jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;        /* 'phi' number in the ring, in {1,4*nr} */
  if(jp > nl4)
    jp = jp - nl4;
  if(jp < 1)
    jp = jp + nl4;

  *phi = (jp - (kshift+1)*0.5) * (piover2 / nr);
}

/* Gives theta and phi corresponding to pixel ipix (RING) for a parameter nside. */
void pix2ang_ring(long nside, long ipix, double *theta, double *phi) {
  int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
  double fact1, fact2, fodd, hip, fihip;
  
  if(nside < 1 || nside > MAX_NSIDE) {
    eprintf("nside out of range: %ld\n", nside);
    exit(0);
  }
  npix = 12*nside*nside;        /* total number of points */
  if(ipix < 0 || ipix > npix-1) {
    eprintf("ipix out of range: %ld\n", ipix);
    exit(0);
  }
  
  ipix1 = ipix + 1; /* in {1, npix} */
  nl2 = 2*nside;
  nl4 = 4*nside;
  ncap = 2*nside*(nside-1);     /* points in each polar cap, =0 for nside =1 */
  fact1 = 1.5*nside;
  fact2 = 3.0*nside*nside;
  
  if(ipix1 <= ncap) {                           /* North Polar cap */
    hip = ipix1/2.;
    fihip = floor(hip);
    iring = (int)floor(sqrt(hip - sqrt(fihip))) + 1;   /* counted from North pole */
    iphi = ipix1 - 2*iring*(iring - 1);
    
    *theta = acos(1. - iring*iring / fact2);
    *phi = (iphi - 0.5) * pi/(2.*iring);
  }
  else if(ipix1 <= nl2*(5*nside+1)) {           /* Equatorial region */
    ip = ipix1 - ncap - 1;
    iring = (int)floor(ip / nl4) + nside;       /* counted from North pole */
    iphi = (ip % nl4) + 1;
    
    fodd = 0.5 * (1 + ((iring+nside) % 2));     /* 1 if iring+nside is odd, 1/2 otherwise */
    *theta = acos((nl2 - iring) / fact1);
    *phi = (iphi - fodd) * pi/(2.*nside);
  }
  else {                                        /* South Polar cap */
    ip = npix - ipix1 + 1;
    hip = ip/2.;
    fihip = floor(hip);
    iring = (int)floor(sqrt(hip - sqrt(fihip))) + 1;    /* counted from South pole */
    iphi = 4*iring + 1 - (ip - 2*iring*(iring-1));
    
    *theta = acos(-1. + iring*iring / fact2);
    *phi = (iphi - 0.5) * pi/(2.*iring);
  }
}

long nside2npix(long nside) {
  return 12*nside*nside;
}

long npix2nside(long npix) {
  return (long)floor(sqrt(npix/12.)+0.5);
}

void ang2vec(double theta, double phi, double *vec) {
  double sz;
  if(theta < 0. || theta > pi) {
    eprintf("theta out of range: %f\n", theta);
    exit(0);
  }
  sz = sin(theta);
  vec[0] = sz * cos(phi);
  vec[1] = sz * sin(phi);
  vec[2] = cos(theta);
}

void vec2ang(double *vec, double *theta, double *phi) {
  double norm = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  *theta = acos(vec[2]/norm);
  *phi = 0.0;
  if(vec[0] != 0.0 || vec[1] != 0.0) {
    *phi = atan2(vec[1],vec[0]);        /* in (-pi, pi] */
    if(*phi < 0.0) *phi += twopi;       /* in  [0, 2pi) */
  }
}

float *read_healpix_map(const char *infile, long *nside, char *coordsys, char *ordering) {
  long naxes, *naxis, npix, npercol, irow;
  int status, hdutype, nfound, anynul;
  float nulval, *map;
  char comment[FLEN_COMMENT];
  fitsfile *fptr;

  status = 0;
  map = NULL;
  naxis = NULL;

  /* Open the file */
  if(fits_open_file(&fptr, infile, READONLY, &status) != 0)
    goto error;

  /* Move to the HDU */
  if(fits_movabs_hdu(fptr, 2, &hdutype, &status) != 0)
    goto error;
  if(hdutype != BINARY_TBL) {
    eprintf("Extension is not binary!\n");
    goto error;
  }

  /* Read the sizes of the array */
  if(fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status) != 0)
    goto error;
  naxis = (long *)malloc(naxes*sizeof(long));
  if(fits_read_keys_lng(fptr, "NAXIS", 1, naxes, naxis, &nfound, &status) != 0 || nfound != naxes)
    goto error;
  if(fits_read_key_lng(fptr, "NSIDE", nside, comment, &status) != 0)
    goto error;
  npix = 12*(*nside)*(*nside);
  if((npix % naxis[1]) != 0) {
    eprintf("Problem with npix.\n");
    goto error;
  }
  npercol = npix/naxis[1];

  if(fits_read_key(fptr, TSTRING, "COORDSYS", coordsys, comment, &status) != 0) {
    fprintf(stderr, "WARNING: Could not find %s keyword in file %s\n", "COORDSYS", infile);
    status = 0;
  }

  if(fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status) != 0) {
    fprintf(stderr, "WARNING: Could not find %s keyword in file %s\n", "ORDERING", infile);
    status = 0;
  }

  /* Read the array */
  map = (float *)malloc(npix*sizeof(float));
  nulval = HEALPIX_NULLVAL;
  for(irow = 0; irow < naxis[1]; irow++) {
    if(fits_read_col(fptr, TFLOAT, 1, irow+1, 1, npercol, &nulval, &map[irow*npercol], &anynul, &status))
      goto error;
  }

  if(fits_close_file(fptr, &status) != 0)
    goto error;

  /* Normal return */
  free(naxis);
  return map;

error:
  /* Proper cleanup on error */
  free(map);
  free(naxis);
  fits_report_error(stderr, status);
  return NULL;
}

int write_healpix_map(const float *signal, long nside, const char *filename, char nest, const char *coordsys) {
  fitsfile *fptr;
  int status, hdutype;
  
  int bitpix = SHORT_IMG;
  long naxis = 0;
  long naxes[] = {0,0};
  
  int tfields = 1;
  long nrows;
  
  char extname[] = "BINTABLE";  /* extension name */
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1E" };
  char *tunit[] = { " " };
  const char *order;            /* HEALPix ordering */
  
  /* Standardize coordinate system string */
  if(coordsys[0] == 'G')
    coordsys = "G       ";
  else if(coordsys[0] == 'E')
    coordsys = "E       ";
  else if(coordsys[0] == 'C' || coordsys[0] == 'Q')
    coordsys = "C       ";
  else {
    eprintf("System Coordinates is not correct (Galactic,Ecliptic,Celestial=Equatorial). Celestial system was set.\n");
    coordsys = "C       ";
  }

  /* Choose ordering */
  if(nest)
    order = "NESTED  ";
  else
    order = "RING    ";

  /* Calculate the number of pixels in the full map */
  nrows = 12*nside*nside;
  
  /* Initialize status before calling fitsio routines */
  status = 0;
  
  /* Create new FITS file */
  if(fits_create_file(&fptr, filename, &status) != 0)
    eprintf("Could not create new fits file.\n");
  
  if(fits_create_img(fptr, bitpix, naxis, naxes, &status) != 0)
    eprintf("Could not create new image file.\n");
  
  if(fits_write_date(fptr, &status) != 0)
    eprintf("Could not add date.\n");
  
  /* Move to first HDU  */
  if(fits_movabs_hdu(fptr, 1, &hdutype, &status)) 
    eprintf("Could not move to first HDU.\n");
  
  /* Append a new empty binary table onto the FITS file */
  if(fits_create_tbl(fptr, BINARY_TBL, nrows, tfields, ttype, tform, tunit, extname, &status))
    eprintf("Could not create new binary table.\n");
  
  if(fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX", "HEALPIX Pixelisation", &status))
    eprintf("Could not write PIXTYPE keyword.\n");
  
  if(fits_write_key(fptr, TSTRING, "ORDERING", (void*)order, "Pixel ordering scheme, either RING or NESTED", &status))
    eprintf("Could not write ORDERING keyword.\n");
  
  if(fits_write_key(fptr, TLONG, "NSIDE", &nside, "Resolution parameter for HEALPIX", &status))  
    eprintf("Could not write NSIDE keyword.\n");
 
  if(fits_write_key(fptr, TSTRING, "COORDSYS", (void*)coordsys, "Pixelisation coordinate system", &status))
    eprintf("Could not write COORDSYS keyword.\n");

  if(fits_write_comment(fptr,"           G = Galactic, E = ecliptic, C = celestial = equatorial  ", &status))
    eprintf("Could not write COORDSYS explanation keyword.\n");
  
  if(fits_write_col(fptr, TFLOAT, 1, 1, 1, nrows, (void*)signal, &status))
    eprintf("Could not write signal.\n");

  if(fits_close_file(fptr, &status))
    eprintf("Could not close file.\n");
  
  return 0;
}
