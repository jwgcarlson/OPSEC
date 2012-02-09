#!/usr/bin/python
#
# binpk.py
#   Re-bin an input P(k) according to a given binning scheme.

import getopt, numpy, sys
from math import *

def usage(status):
    f = (status == 0) and sys.stdout or sys.stderr
    print >> f, "Usage: %s [-b bins] [-w m] infile [outfile]" % sys.argv[0]
    print >> f, "Binning:"
    print >> f, "  (Not all of this is implemented yet.)"
    print >> f, "  k0,k1,...,kN     -- N bins given by intervals [k0,k1], [k1,k2], etc."
    print >> f, "  kmin-kmax:N      -- N linearly spaced bins between kmin and kmax"
    print >> f, "  kmin-kmax:N:log  -- N logarithmically spaced bins (kmin must be > 0)"
    print >> f, "Weighting:"
    print >> f, "  The argument to the -w switch is a real number 'm'.  The program computes"
    print >> f, "  the average value of k^m P(k) within each bin.  For binning a 3D power"
    print >> f, "  spectrum, choose m = 2 (although for narrow bins this makes little difference."
    print >> f, "Output file defaults to stdout."
    sys.exit(status)

if __name__ == '__main__':
    verbose = False
    bins = []
    m = 2.
    numfmt = "%e"

    # Parse command line options
    (opts,args) = getopt.getopt(sys.argv[1:], "hvb:f:w:")
    for (o,a) in opts:
        if o == '-h':
            usage(0)
        elif o == '-v':
            verbose = True
        elif o == '-b':
            # The argument to the '-b' option is a comma-separated list of bin
            # specification.  A bin specification is of the form
            #    k1-k2   or   k1-k2:N
            # where k1 and k2 are positive real values, and N is a positive
            # integer.  This specifies N equally spaced bins between k1 and
            # k2.  In the first form, N defaults to 1.
            for bs in a.split(','):
                try:
                    fields = bs.split(':')
                    kk = fields[0].split('-')
                    k1, k2 = float(kk[0]), float(kk[1])
                    N = (len(fields) >= 2) and int(fields[1]) or 1
                    mapping = (len(fields) >= 3) and fields[2] or "linear"
                    if mapping.startswith("lin"):
                        identity = lambda x: x
                        f, finv = identity, identity
                    elif mapping.startswith("log"):
                        f, finv = log, exp
                    fk1, fk2 = f(k1), f(k2)
                    for i in range(N):
                        b = ( finv(fk1 + i*(fk2 - fk1)/N), finv(fk1 + (i+1)*(fk2 - fk1)/N) )
                        bins.append(b)
                except Exception, e:
                    print >> sys.stderr, "Error parsing binning specification '%s'" % bs
                    usage(1)
        elif o == '-f':
            pass
        elif o == '-w':
            pass
        else:
            print >> sys.stderr, "Unrecognized option: %s" % o
            usage(1)

    # Make sure we have at least an input file on the command line (output file defaults to stdout)
    if len(args) < 1:
        print >> sys.stderr, "Must specify input P(k) file"
        usage(1)
    elif len(args) > 3:
        print >> sys.stderr, "Too many arguments provided"
        usage(1)

    # Make sure we can read from the input file
    fin = open(args[0], 'r')
    if not fin:
        print >> sys.stderr, "Could not open %s for reading" % fin.name
        sys.exit(1)
    if verbose:
        print "Reading input P(k) from %s" % fin.name

    # Make sure we can write to the output file
    fout = (len(args) < 2) and sys.stdout or open(args[1], 'w')
    if not fout:
        print >> sys.stderr, "Could not open %s for writing" % fout.name
        sys.exit(1)
    if verbose:
        print "Writing binned P(k) to %s" % fout.name

    # Make sure we have reasonable bins
    if bins == []:
        print >> sys.stderr, "Bins not specified, assuming 20 linear bins over [0, 0.2]"
        bins = [(k,k+0.01) for k in numpy.arange(0.0, 0.2, 0.01)]
    else:
#        print "(debugging) bins = %s" % bins
        pass

    # Read in data
    data = numpy.loadtxt(fin)
    kk = data[:,0]
    nk = len(kk)
    pk = data[:,1:]
    npk = pk.shape[1]

    nbins = len(bins)
    binned_pk = numpy.zeros((nbins,npk))
    counts = numpy.zeros(nbins)

    # Accumulate P(k) in bins
    for i in range(nk):
        k = kk[i]
        for j in range(nbins):
            (kmin,kmax) = bins[j]
            if kmin <= k < kmax:
                binned_pk[j,:] += k**m * pk[i,:]
                counts[j] += k**m

    # Compute average value within bin
    for j in range(nbins):
        binned_pk[j,:] /= counts[j]

    # Write results to file
    for j in range(nbins):
        (kmin,kmax) = bins[j]
        print >> fout, numfmt % kmin,
        for c in range(npk):
            print >> fout, (" " + numfmt) % binned_pk[j,c],
        print >> fout, ""
        print >> fout, numfmt % kmax,
        for c in range(npk):
            print >> fout, (" " + numfmt) % binned_pk[j,c],
        print >> fout, ""
