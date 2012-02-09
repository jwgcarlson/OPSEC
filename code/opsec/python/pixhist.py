#!/usr/bin/python

import numpy
import matplotlib.mlab as mlab
import matplotlib.pyplot as pyplot

from math import *
from abn import *
from cfg import *

cfg = cfg_read("klpk.cfg")
Nmodes = int(cfg["Nmodes"])
evalfile = cfg["evalfile"]
pixfile = cfg["valfile"]

evals = numpy.array(abn_read(evalfile))
y = numpy.array(abn_read(pixfile))
C = 1 + evals   # diagonal of covariance matrix
z = y/numpy.sqrt(C)

n, bins, patches = pyplot.hist(z, 50, normed=1, facecolor='green', alpha=0.75)

# Best fit line: unit Gaussian
z = mlab.normpdf(bins, 0., 1.)
pyplot.plot(bins, z, 'r--', linewidth=1)

pyplot.xlabel('y_i/(C_{ii})^{1/2}')
pyplot.ylabel('Distribution')
pyplot.axis([-3, 3, 0, 1])

pyplot.show()
