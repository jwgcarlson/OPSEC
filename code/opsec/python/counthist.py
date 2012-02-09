#!/usr/bin/python

import numpy, sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot

from math import *
from abn import *
from cfg import *

cfg = cfg_read("klpk.cfg")
countfile = cfg["countfile"]

data = numpy.loadtxt(countfile)
N = data[:,1]
Nbar = data[:,2]
var = data[:,3]
x = (N - Nbar)/numpy.sqrt(var)

n, bins, patches = pyplot.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

# Best fit line: unit Gaussian
z = mlab.normpdf(bins, 0., 1.)
pyplot.plot(bins, z, 'r--', linewidth=1)

pyplot.xlabel('x_a/(1 + S_{aa})^{1/2}')
pyplot.ylabel('Distribution')
pyplot.axis([-3, 3, 0, 1])

pyplot.show()
