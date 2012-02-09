#!/usr/bin/python

import sys
from numpy import array, concatenate, diag, dot, identity, linspace, loadtxt, ones, sqrt, trace, transpose, zeros
from numpy.linalg import eig, inv, svd

from abn import *
from cfg import *


cfg = cfg_read("klpk.cfg")
opts = {}

# Read in signal matrix
data = array(abn_read(cfg["evalfile"], opts))
Nmodes = int(opts["Nmodes"])
S = diag(data)
N = identity(Nmodes)
C = S + N
Cinv = inv(C)

# Read in parameter derivatives of the covariance matrix
data = array(abn_read(cfg["covfile"], opts))
Nparams = int(opts["Nparams"])
assert Nmodes == int(opts["Nmodes"])
Cp = data.reshape((Nparams,Nmodes,Nmodes))

# Read in pixel values, y_i
y = array(abn_read(cfg["valfile"], opts))

# Compute Fisher matrix
F = zeros((Nparams,Nparams), dtype=float)
for m in range(Nparams):
    for n in range(Nparams):
        F[m,n] = 0.5 * trace(dot(Cinv, dot(Cp[m], dot(Cinv, Cp[n]))))

def sqrtm(M):
    """Compute square root of positive-definite symmetric matrix."""
    (U,s,Vh) = svd(M)
    D = diag(sqrt(s))
    return dot(U, dot(D, Vh))

# Choose M matrix
#M = inv(F)
#M = identity(Nparams)
M = inv(sqrtm(F))

# ...and normalize it
W = dot(M,F)
for m in range(Nparams):
    M[m,:] /= sum(W[m,:])
W = dot(M,F)
#print M

# Compute quadratic estimators
q = zeros(Nparams)
for m in range(Nparams):
    q[m] = 0.5 * dot(y, dot(Cinv, dot(Cp[m], dot(Cinv, y))))
#print q

# Compute power estimates
phat = dot(M, q)

# Compute error bars
variance = dot(M, dot(F, transpose(M)))
sigma = sqrt(diag(variance))


try:
    # Try to plot the estimate and the prior
    from matplotlib import pyplot

    # Read in model parameters
    pkfile = cfg["RealBandedModel.pkfile"]
    Nbands = int(cfg["RealBandedModel.Nbands"])
    kmin = float(cfg["RealBandedModel.kmin"])
    kmax = float(cfg["RealBandedModel.kmax"])

    # Load prior P(k)
    data = loadtxt(pkfile)
    kpri,ppri = data[:,0], data[:,1]

    # Prepare data for plotting
    kbands = linspace(kmin, kmax, Nbands+1)
    left = kbands[:Nbands]
    width = (kmax - kmin)/Nbands

    # Plot estimate over true P(k)
    pyplot.bar(left, phat, width, yerr=sigma)
    pyplot.plot(kpri, ppri, 'r--')
    pyplot.xlim(0, 0.25)
    pyplot.show()
except:
    print "Couldn't import matplotlib."
    print "phat = %s" % phat
