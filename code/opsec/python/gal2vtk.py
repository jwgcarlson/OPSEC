#!/usr/bin/python

import os.path
import struct
import sys
from collections import namedtuple
from math import *

from abn import *
from cfg import *


def convert_galaxies_to_vtk(galfile, outfile=None, binary=True):
    binary = False      # binary output currently causes Paraview to segfault, need to debug

    # Read in galaxies
    opts = {}
    galaxies = abn_read(galfile, opts)
    if not galaxies:
        print >> sys.stderr, "Error reading in galaxies from '%s'" % galfile
        return 1
    Ngals = len(galaxies)

    # Decide whether we're in spherical or Cartesian coordinates
    if not "coordsys" in opts:
        print "Coordinate system not specified, defaulting to spherical"
        coordmap = lambda (r,mu,phi): (r*sqrt(1-mu**2)*cos(phi), r*sqrt(1-mu**2)*sin(phi), r*mu)
    if opts["coordsys"] == "spherical":
        coordmap = lambda (r,mu,phi): (r*sqrt(1-mu**2)*cos(phi), r*sqrt(1-mu**2)*sin(phi), r*mu)
    elif opts["coordsys"] == "cartesian":
        coordmap = lambda (x,y,z): (x,y,z)
    else:
        print >> sys.stderr, "Error, unrecognized value for 'coordsys': %s" % opts["coordsys"]
        return 1

    # Open output file
    if outfile is None:
        # Use the same basename as galfile, but with a .vtk extension
        root, ext = os.path.splitext(galfile)
        outfile = root + ".vtk"
    fp = open(outfile, "wb")
    if not fp:
        print >> sys.stderr, "Error, could not open output file '%s' for writing" % outfile

    # Write VTK header
    print >> fp, "# vtk DataFile Version 2.0"
    print >> fp, "Visualization of galaxies from '%s'" % galfile
    print >> fp, binary and "BINARY" or "ASCII"
    print >> fp, "DATASET UNSTRUCTURED_GRID"

    # Write points
    print >> fp, "POINTS %d float" % Ngals
    for p in galaxies:
        x = coordmap(p)
        if binary:
            fp.write(struct.pack("=3f", x[0], x[1], x[2]))
        else:
            print >> fp, "%f %f %f" % x

    # Write cells
    print >> fp, "CELLS %d %d" % (Ngals,2*Ngals)
    for i in range(Ngals):
        if binary:
            fp.write(struct.pack("=2i", 1, i))
        else:
            print >> fp, "1 %d" % i

    # Write cell types
    VTK_VERTEX = 1
    print >> fp, "CELL_TYPES %d" % Ngals
    for i in range(Ngals):
        if binary:
            fp.write(struct.pack("=i", VTK_VERTEX))
        else:
            print >> fp, VTK_VERTEX

    fp.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print >> sys.stderr, "Usage: %s <galfile>" % (sys.argv[0])
        sys.exit(1)

    galfile = sys.argv[1]
    convert_galaxies_to_vtk(galfile, binary=False)
