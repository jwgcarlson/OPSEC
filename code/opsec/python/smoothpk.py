#!/usr/bin/python
#
# smoothpk.py:
#   Smooth a given linear P(k) to obtain a desired variance-in-cells.

import getopt, os, sys
import numpy as np
from math import *
from abn import *
from cfg import *

def usage():
    print "Usage: %s [-c cfgfile]" % sys.argv[0]

if __name__ == "__main__":
    cfg = {}

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
            cfg.update(cfg_read(a))
        else:
            assert False, "unhandled option '%s' (this should not happen)" % o

    # Parse configuration options
    for a in args:
        cfg.update(cfg_read_line(a))

    if "infile" not in cfg:
        print "Missing configuration option 'infile'"
        sys.exit(1)
    infile = cfg["infile"]

    kcol, pcol = 0, 1
    if "kcol" in cfg: kcol = int(cfg["kcol"]) - 1
    if "pcol" in cfg: pcol = int(cfg["pcol"]) - 1

    if "outfile" not in cfg:
        print "Missing configuration option 'outfile'"
        sys.exit(1)
    outfile = cfg["outfile"]

    if "N" not in cfg:
        print "Missing configuration option 'N'"
        sys.exit(1)
    N = float(cfg["N"])

    if "L" not in cfg:
        print "Missing configuration option 'L'"
        sys.exit(1)
    L = float(cfg["L"])

    p = 0.01
    if "p" in cfg: p = float(cfg["p"])

    # Find sigma for which Pr(\delta < -1) = p
    def Phi(x): return 0.5*erfc(-x/sqrt(2))
    sigma = 0.01
    while Phi(-1/sigma) < p:
        sigma += 0.0001
    print "desired sigma = %g" % sigma

    # Load power spectrum
    data = np.loadtxt(infile)
    kdata = data[:,kcol]
    pdata = data[:,pcol]
    def P(k): return np.interp(k, kdata, pdata)

    knyq = pi*N/L
    Nk = 8192
    k = np.linspace(0, knyq, Nk+1)
    k2 = k**2
    A = k2*P(k) / (2*pi**2)

    def ComputeSigma(R):
        return sqrt( np.trapz(A * np.exp(-k2*R**2), dx=knyq/Nk) )

    # Find smoothing radius R for which <\delta^2> = \sigma^2
    R = 0.0
    while ComputeSigma(R) > sigma:
        R += 0.01
    print "R = %g" % R

    # Modify data and save to file
    pdata *= np.exp(-kdata**2*R**2)
    np.savetxt(outfile, data, fmt='%e')
