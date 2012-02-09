#!/usr/bin/python
#
import sys
from collections import namedtuple
from math import *
from pyx import *

from abn import *
from cfg import *

cfgfile = "opsec.cfg"

def read_cells(cellfile):
    Cell = namedtuple("Cell", "a G rmin rmax mumin mumax phimin phimax Veff Nbar")
    cells = abn_read(cellfile)
    return [Cell(*c) for c in cells]

def read_galaxies(galfile):
    Galaxy = namedtuple("Galaxy", "r mu phi")
    gals = abn_read(galfile)
    return [Galaxy(*g) for g in gals]

def hammer(mu, phi):
    """Return the Hammer-Aitoff projection (x,y) of the coordinates (mu,phi) on the sphere."""
    mu = min(+1, max(-1, mu))   # clamp mu to [-1,+1]
#    phi += 2*pi*((phi < 0) - (phi >= 2*pi))      # choose phi in [0,2*pi)
    sintheta = sqrt(1 - mu**2)
    lon = phi - pi
    z = sqrt(1 + sintheta * cos(lon/2))
    x = 2*sintheta * sin(lon/2) / z
    y = mu / z
    return (5*x,5*y)

def hammer_path(mu0, phi0, mu3, phi3):
    """Create a curved path from the current position (mu0,phi0) to (mu3,phi3)."""
    mu1 = (2/3.)*mu0 + (1/3.)*mu3
    mu2 = (1/3.)*mu0 + (2/3.)*mu3
    phi1 = (2/3.)*phi0 + (1/3.)*phi3
    phi2 = (1/3.)*phi0 + (2/3.)*phi3
    x1,y1 = hammer(mu1,phi1)
    x2,y2 = hammer(mu2,phi2)
    x3,y3 = hammer(mu3,phi3)
    return path.curveto(x1,y1, x2,y2, x3,y3)

def hammer_outline():
    """Create a path outlining the Hammer-Aitoff projection shape."""
    r = path.path(path.moveto(*hammer(1,2*pi)))
    mus = [1 - 2*i/100. for i in range(1,100)]
    for mu in mus:
        r.append(path.lineto(*hammer(mu, 2*pi)))
    r.append(path.lineto(*hammer(-1,2*pi)))
    mus.reverse()
    for mu in mus:
        r.append(path.lineto(*hammer(mu, 0)))
    r.append(path.closepath())
    return r

def draw_cell(c, p, style):
    x0,y0 = hammer(p.mumin, p.phimin)
    r = path.path(path.moveto(x0,y0))
    r.append(hammer_path(p.mumin, p.phimin, p.mumax, p.phimin))
    r.append(hammer_path(p.mumax, p.phimin, p.mumax, p.phimax))
    r.append(hammer_path(p.mumax, p.phimax, p.mumin, p.phimax))
    r.append(path.closepath())
    c.stroke(r, style)


if __name__ == '__main__':
    cfg = cfg_read(cfgfile)
    coordsys = cfg["coordsys"]
    if coordsys != "spherical":
        print >> sys.stderr, "gals2pdf.py: only makes sense for spherical coordinates"
        sys.exit(1)
    cellfile = cfg["cellfile"]
    galfile = cfg["galfile"]
    Nr = int(cfg["Nr"])
    RMin = float(cfg["RMin"])
    RMax = float(cfg["RMax"])
    MuMin = float(cfg["MuMin"])
    MuMax = float(cfg["MuMax"])
    PhiMin = float(cfg["PhiMin"])
    PhiMax = float(cfg["PhiMax"])

    cells = read_cells(cellfile)
    gals = read_galaxies(galfile)

    pages = []
    for d in range(Nr):
        rmin = RMin + d*(RMax - RMin)/Nr
        rmax = RMin + (d+1)*(RMax - RMin)/Nr

        c = canvas.canvas()
        c.stroke(hammer_outline(), [style.linewidth.normal, color.rgb.black])

        # Draw cell outlines in gray
        for p in cells:
            if rmin < 0.5*(p.rmin + p.rmax) < rmax:
                draw_cell(c, p, [style.linewidth.Thin, color.gradient.Gray.getcolor(0.25)])

        # Draw galaxies, with color representing depth (blue for near, red for far)
        # TODO: figure out how to make dots always appear as 1 pixel
        for g in gals:
            if rmin <= g.r < rmax:
                mu = g.mu
                phi = g.phi + 2*pi*((g.phi < 0) - (g.phi >= 2*pi))      # choose phi in [0,2*pi)
                x, y = hammer(mu, phi)
                u = (g.r - RMin)/(RMax - RMin)
                c.stroke(path.line(x, y, x+0.04, y), [style.linewidth.Thin, color.gradient.BlueRed.getcolor(u)])

        pages.append(document.page(c))

    doc = document.document(pages)
    doc.writePDFfile("gals.pdf")
