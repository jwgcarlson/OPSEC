#!/usr/bin/python
#
# plot_estimate.py:
#   Plot the estimated P(k).

import getopt, os, sys
import numpy as np
from abn import *
from cfg import *
from matplotlib import pylab

# Default configuration location
cfgfile = "opsec.cfg"


def usage():
    print "Usage: %s [-c config]" % sys.argv[0]

def read_estimates(estfile):
    # Open parameter estimates file
    fest = open(estfile, "rb")
    if not fest:
        print "read_estimates: could not open '%s'" % estfile
        sys.exit(1)

    # Read estimates and related matrices
    opts = {}
    phat = abn_read_numpy(fest, opts)
    Nparams = int(opts["Nparams"])
    cov = abn_read_numpy(fest).reshape((Nparams,Nparams))
    F = abn_read_numpy(fest).reshape((Nparams,Nparams))
    W = abn_read_numpy(fest).reshape((Nparams,Nparams))
    M = abn_read_numpy(fest).reshape((Nparams,Nparams))
    q = abn_read_numpy(fest)
    f = abn_read_numpy(fest)

    fest.close()
    return (Nparams,phat,cov,F,W,M,q,f)

def plot_real_model(plotfile, cfg, phat, cov, W):
    # Read in model parameters
    pkfile = cfg["RealModel.pkfile"]
    Nbands = int(cfg["RealModel.Nbands"])
    kmin = float(cfg["RealModel.kmin"])
    kmax = float(cfg["RealModel.kmax"])
    assert len(phat) == Nbands

    # Plot y = 0 as a baseline
    k = np.array([0., 1.1*kmax])
    p = np.array([0., 0.])
    pylab.plot(k, p, 'k:')

    # Plot estimate over true P(k)
    k = np.array([kmin + (i+0.5)*(kmax-kmin)/Nbands for i in range(Nbands)])
    p = phat
    sigma = np.sqrt(np.diag(cov))
    pylab.errorbar(k, p, yerr=sigma, fmt='k+', label='OPSEC')

    # Plot prior P(k) in red
    data = np.loadtxt(pkfile)
    k,p = data[:,0], data[:,1]
    pylab.plot(k, p, 'r--', label='prior')

    # Use prior to determine plot bounds
    pmax = p.max()

    # Plot true P(k) used to generate mock catalog as green line (if available)
    if "pktrue" in cfg:
        pktrue = cfg["pktrue"]
        if os.path.exists(pktrue):
            data = np.loadtxt(pktrue)
            pylab.plot(data[:,0], data[:,1], 'g-', label='true')


    # Plot exact P(k) band powers as black lines (if available)
    if "pkexact" in cfg:
        pkexact = cfg["pkexact"]
        if os.path.exists(pkexact):
            data = np.loadtxt(pkexact)
            pylab.plot(data[:,0], data[:,1], 'b-', label='exact')

    pylab.xlim(0, kmax + 0.01)
    pylab.ylim(-0.1*pmax, 1.1*pmax)
    pylab.xlabel('k [h/Mpc]')
    pylab.ylabel('P(k) [(Mpc/h)^3]')
    pylab.legend(loc='upper right')
    pylab.savefig(plotfile)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:", ["help", "config="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)

    # Parse command line switches
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-c", "--config"):
            cfgfile = a
        else:
            assert False, "unhandled option '%s' (this should not happen)" % o

    # Parse configuration options
    if not os.path.exists(cfgfile):
        print "Could not find config file '%s'" % cfgfile
        sys.exit(1)
    cfg = cfg_read(cfgfile)
    for a in args:
        cfg.update(cfg_read_line(a))

    if "estfile" not in cfg:
        print "Missing 'estfile' option in '%s'" % cfgfile
        sys.exit(1)
    estfile = cfg["estfile"]

    if "model" not in cfg:
        print "Missing 'model' option in '%s'" % cfgfile
        sys.exit(1)
    model = cfg["model"]

    if "mixing" in cfg:
        mixing = cfg["mixing"]
    else:
        mixing = "inverse"

    if "plotfile" in cfg:
        plotfile = cfg["plotfile"]
    else:
        plotfile = "pk.pdf"

    (Nparams,phat,cov,F,W,M,q,f) = read_estimates(estfile)

    if model == "RealModel":
        plot_real_model(plotfile, cfg, phat, cov, W)
