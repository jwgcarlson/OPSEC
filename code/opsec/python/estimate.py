#!/usr/bin/python

import getopt, sys
from numpy import array, concatenate, diag, dot, identity, linspace, loadtxt, ones, sqrt, sum, trace, transpose, zeros
from numpy import linalg
from abn import *
from cfg import *

import time

last = time.time()
def dtime():
    global last
    now = time.time()
    dt = now - last
    last = now
    return dt

if __name__ != "__main__":
    raise RuntimeError, "This is a script, not a module."

def usage():
    print "Usage: Use the source, Luke."

try:
    opts, args = getopt.getopt(sys.argv[1:], "hc:", ["help", "config="])
except getopt.GetOptErr, err:
    print str(err)
    usage()
    sys.exit(1)

# Parse command line switches
cfgfile = None
for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-c", "--config"):
        cfgfile = a
    else:
        assert False, "unhandled option"

# Parse configuration options
if cfgfile:
    cfg = cfg_read(cfgfile)
else:
    cfg = {}
for a in args:
    cfg.update(cfg_read_line(a))

# Check for complete configuration
try:
    evalfile = cfg["evalfile"]
    covfile = cfg["covfile"]
    pixfile = cfg["valfile"]
except KeyError:
    print >> sys.stderr, "Missing configuration options."
    raise

# Read in signal matrix
opts = {}
#evals = array(abn_read(evalfile, opts))
evals = abn_read_numpy(evalfile, opts)
Nmodes = int(opts["Nmodes"])
S = diag(evals)
N = identity(Nmodes)
C = S + N
Cinv = diag(1/(1+evals))

print "Reading signal matrix: dt = %g s" % dtime()

# Read in parameter derivatives of the covariance matrix
opts = {}
#data = array(abn_read(covfile, opts))
data = abn_read_numpy(covfile, opts)
Nparams = int(opts["Nparams"])
assert Nmodes == int(opts["Nmodes"])
Cp = data.reshape((Nparams,Nmodes,Nmodes))

print "Reading covariance matrix derivatives: dt = %g s" % dtime()

# Read in pixel values, y_i
opts = {}
#y = array(abn_read(pixfile, opts))
y = abn_read_numpy(pixfile, opts)

print "Reading pixel values: dt = %g s" % dtime()

# Compute Fisher matrix: F_{mn} = (1/2) Tr[C^{-1} C_{,m} C^{-1} C_{,n}]
D = [dot(Cinv, Cp[m]) for m in range(Nparams)]
Dt = [transpose(D[m]) for m in range(Nparams)]
F = zeros((Nparams,Nparams), dtype=float)
for m in range(Nparams):
    for n in range(m+1):
        F[m,n] = 0.5 * sum(D[m]*Dt[n])   # Tr[D_m D_n] = \sum_{ij} (D_m)_{ij} (D_n)_{ji}
for m in range(Nparams):
    for n in range(m+1,Nparams):
        F[m,n] = F[n,m]

print "Computing Fisher matrix: dt = %g s" % dtime()

def sqrtm(M):
    """Compute square root of positive-definite symmetric matrix."""
    (U,s,Vh) = linalg.svd(M)
    D = diag(sqrt(s))
    return dot(U, dot(D, Vh))

# Choose M matrix
#M = linalg.inv(F)
#M = identity(Nparams)
M = linalg.inv(sqrtm(F))

print "Computing M matrix: dt = %g s" % dtime()

# ...and normalize it
W = dot(M,F)
for m in range(Nparams):
    M[m,:] /= sum(W[m,:])
W = dot(M,F)

print "Normalizing M matrix: dt = %g s" % dtime()

# Compute quadratic estimators: q_n = (1/2) y^T C^{-1} C_{,n} C^{-1} y
q = zeros(Nparams)
for n in range(Nparams):
    q[n] = 0.5 * dot(y, dot(D[n], y/(1+evals))) # (C^{-1} y)_i = y_i/(1 + \lambda_i)

print "Computing quadratic estimators: dt = %g s" % dtime()

# Compute shot noise corrections (implicitly using N = Identity)
f = zeros(Nparams)
for n in range(Nparams):
    f[n] = 0.5 * sum(Cinv*transpose(D[n]))      # Tr[C^{-1} N C^{-1} C_{,n}] = \sum_{ij} (C^{-1})_{ij} (D_n)_{ji}

print "Computing shot noise corrections: dt = %g s" % dtime()

# Compute power estimates
phat = dot(M, q - f)
cov = dot(M, dot(F, transpose(M)))

print "Computing power estimates: dt = %g s" % dtime()

# Print relevant information
print "phat = %s" % repr(phat)
print "cov = %s" % repr(cov)
print "W = %s" % repr(W)
print "F = %s" % repr(F)

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
    sigma = sqrt(diag(cov))
    pyplot.bar(left, phat, width, yerr=sigma)
    pyplot.plot(kpri, ppri, 'r--')
    pyplot.xlim(0, 0.25)
    pyplot.show()
except:
    pass
