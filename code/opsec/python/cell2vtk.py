#!/usr/bin/python

import os.path
import struct
import sys
from collections import namedtuple
from math import *

from abn import *
from cfg import *


def convert_cells_to_vtk(cellfile, modefile, modenum, outfile=None, binary=True):
    binary = False      # binary output currently causes Paraview to segfault, need to debug

    # Read in cells
    opts = {}
    cells = abn_read(cellfile, opts)
    if not cells:
        print >> sys.stderr, "Error reading in cells from '%s'" % cellfile
        return 1
    Ncells = len(cells)

    # Read in desired mode
    mode = abn_read(modefile, opts, first=modenum*Ncells, last=(modenum+1)*Ncells-1)
    if not mode:
        print >> sys.stderr, "Error reading in mode %d from '%s'" % (modenum,modefile)
        return 1

    # Decide whether we're in spherical or Cartesian coordinates
    if not "coordsys" in opts:
        print >> sys.stderr, "Error, cells file does not contain the 'coordsys' option"
        return 1
    if opts["coordsys"] == "spherical":
        coordmap = lambda (r,mu,phi): (r*sqrt(1-mu**2)*cos(phi), r*sqrt(1-mu**2)*sin(phi), r*mu)
    elif opts["coordsys"] == "cartesian":
        coordmap = lambda (x,y,z): (x,y,z)
    else:
        print >> sys.stderr, "Error, unrecognized value for 'coordsys': %s" % opts["coordsys"]
        return 1

    # Open output file
    if outfile is None:
        # Use the same basename as cellfile, but with a .vtk extension
        root, ext = os.path.splitext(cellfile)
        outfile = root + ".vtk"
    fp = open(outfile, "wb")
    if not fp:
        print >> sys.stderr, "Error, could not open output file '%s' for writing" % outfile

    points = []         # list of all points in unstructured grid (vertices of cells)
    pointmap = {}       # map between triple (x,y,z) and index i into 'pts' array
    vtkcells = []

    Cell = namedtuple("Cell", "a G x1 x2 y1 y2 z1 z2 Veff Nbar")
    cells = map(Cell._make, cells)
    VTK_HEXAHEDRON = 12
    for c in cells:
        indices = []
        for vertex in [(c.x1,c.y1,c.z1), (c.x1,c.y1,c.z2), (c.x1,c.y2,c.z2), (c.x1,c.y2,c.z1),
                       (c.x2,c.y1,c.z1), (c.x2,c.y1,c.z2), (c.x2,c.y2,c.z2), (c.x2,c.y2,c.z1)]:
            p = coordmap(vertex)
            if p in pointmap:
                i = pointmap[p]
            else:
                i = len(points)
                points.append(p)
                pointmap[p] = i
            indices.append(i)
        vtkcells.append((indices, VTK_HEXAHEDRON))

    # Write VTK header
    print >> fp, "# vtk DataFile Version 2.0"
    print >> fp, "Visualization of cells from '%s'" % cellfile
    print >> fp, binary and "BINARY" or "ASCII"
    print >> fp, "DATASET UNSTRUCTURED_GRID"

    # Write points
    print >> fp, "POINTS %d float" % len(points)
    for p in points:
        if binary:
            fp.write(struct.pack("=3f", p[0], p[1], p[2]))
        else:
            print >> fp, "%f %f %f" % p

    # Write cells
    size = sum([1 + len(vc[0]) for vc in vtkcells])
    print >> fp, "CELLS %d %d" % (Ncells,size)
    for vc in vtkcells:
        indices = vc[0]
        n = len(indices)
        if binary:
            fp.write(struct.pack("=i", n))
            fp.write(struct.pack("=%di" % n, *indices))
        else:
            print >> fp, n,
            for i in indices:
                print >> fp, i,
            print >> fp, ""

    # Write cell types
    print >> fp, "CELL_TYPES %d" % Ncells
    for vc in vtkcells:
        ct = vc[1]
        if binary:
            fp.write(struct.pack("=i", ct))
        else:
            print >> fp, ct

    # Write mode values as a scalar attribute
    print >> fp, "CELL_DATA %d" % Ncells
    print >> fp, "SCALARS psi float 1"
    print >> fp, "LOOKUP_TABLE default"
    for val in mode:
        if binary:
            fp.write(struct.pack("=f", val))
        else:
            print >> fp, val

    fp.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print >> sys.stderr, "Usage: %s <cellfile> <modefile> <modenum>" % (sys.argv[0])
        sys.exit(1)

    cellfile = sys.argv[1]
    modefile = sys.argv[2]
    modenum = int(sys.argv[3])
    convert_cells_to_vtk(cellfile, modefile, modenum, binary=False)
