#!/usr/bin/python
#
import sys
from collections import namedtuple
from math import *
from pyx import *

from abn import *
from cfg import *


def read_modes(modefile):
    opts = {}
    modes = abn_read_numpy(modefile, opts)
    n = int(opts["Ncells"])
    nev = int(opts["Nmodes"])
    return modes.reshape((nev,n))

def read_cells(cellfile):
    Cell = namedtuple("Cell", "a G rmin rmax mumin mumax phimin phimax Veff Nbar")
    cells = abn_read(cellfile)
    return [Cell(*c) for c in cells]


def hammer(mu, phi):
    """Return the Hammer-Aitoff projection (x,y) of the coordinates (mu,phi) on the sphere."""
    if mu > 1: mu = 1.
    elif mu < -1: mu = -1.
    sintheta = sqrt(1 - mu**2)
    lon = phi - pi
    z = sqrt(1 + sintheta * cos(lon/2))
    x = 2*sintheta * sin(lon/2) / z
    y = mu / z
    return (5*x,5*y)

def hammer_path(mu0, phi0, mu3, phi3):
    mu1 = (2/3.)*mu0 + (1/3.)*mu3
    mu2 = (1/3.)*mu0 + (2/3.)*mu3
    phi1 = (2/3.)*phi0 + (1/3.)*phi3
    phi2 = (1/3.)*phi0 + (2/3.)*phi3
    x1,y1 = hammer(mu1,phi1)
    x2,y2 = hammer(mu2,phi2)
    x3,y3 = hammer(mu3,phi3)
    return path.curveto(x1,y1, x2,y2, x3,y3)

def hammer_outline():
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
    c.fill(r, style)

#def draw_angular_mode(outfile, cells, values):
def draw_angular_mode(cells, values):
    c = canvas.canvas()
    c.stroke(hammer_outline())

    n = len(values)
    v = [cells[i].Nbar**0.5/cells[i].Veff * values[i] for i in range(n)]
    vmax = max([abs(v[i]) for i in range(n)])
    red_gradient = color.lineargradient(color.rgb.white, color.rgb.red)
    blue_gradient = color.lineargradient(color.rgb.white, color.rgb.blue)
    print "vmax = %f" % (vmax)
    for i in range(n):
        if v[i] >= 0:
            draw_cell(c, cells[i], [red_gradient.getcolor(v[i]/vmax)])
        else:
            draw_cell(c, cells[i], [blue_gradient.getcolor(-v[i]/vmax)])

#    c.writeEPSfile(outfile)
#    c.writePDFfile(outfile)
    return document.page(c)



if __name__ == '__main__':
    cells = read_cells("cells.abn")
    evals = abn_read_numpy("evals.abn")
    modes = read_modes("modes.abn")

    # Select all cells in the innermost shell
    r = cells[0].rmin

    pages = []
    for i in range(5):
        shell_cells = [p for p in cells if p.rmin == r]
        shell_values = [modes[i][p.a] for p in shell_cells]
        pages.append(draw_angular_mode(shell_cells, shell_values))
    d = document.document(pages)
    d.writePDFfile("angular.pdf")
