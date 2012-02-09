# abn.py
#
# Python routines to read and write .abn files.
# Currently implemented routines:
# - abn_read(f, opts=None, first=0, last=-1)
# - abn_read_numpy(f, opts=None, first=0, last=-1)

# TODO:
# - add documentation
# - implement abn_write()
# - add object-oriented interface
# - allow data blocks to be read a few entries at a time
# - test non-native endianness
# - add abn_read_header()

import re, struct, sys
from cfg import cfg_read_line

SEEK_CUR = 1
begindata_pattern = re.compile(r"#!Begin data block: n = (\d+), size = (\d+), endian = (\w), format = (\w+)")

def abn_read(f, opts=None, first=0, last=-1):
    """Read a data block from a .abn file.  The file pointer is placed after
    the last element read."""

    # If f is a file object, fine.  Otherwise assume it is a filename, and open the file.
    if isinstance(f, file):
        filename = f.name
    else:
        filename = f
        f = open(filename, "rb")
    if not f:
        print >> sys.stderr, "abn_read: could not open '%s'" % filename
        return None

    # Read .abn header
    n = 0
    size = 0
    endian = ''
    fmt = ''
    line = f.readline()
    while line:
        m = re.match(begindata_pattern, line)
        if m:
            n = int(m.group(1))
            size = int(m.group(2))
            endian = m.group(3)
            fmt = m.group(4) 
            break
        elif opts is not None:
            opts.update(cfg_read_line(line))
        line = f.readline()

    if not line:        # reached end of file without finding a data header
        print >> sys.stderr, "abn_read: could not read data header from '%s'" % filename
        return None

    # Check for sensible limits on requested data
    if first < 0 or first > n-1:
        first = 0
    if last < 0  or last > n-1:
        last = n-1

    # Read binary data
    data = []
    assert size == struct.calcsize(fmt)
    if endian == 'L':
        fmt = '<' + fmt
    else:
        fmt = '>' + fmt
    f.seek(first*size, SEEK_CUR)     # skip ahead to first requested entry
    for i in range(first,last+1):
        s = f.read(size)
        entry = struct.unpack(fmt, s)
        if len(entry) == 1:
            data.append(entry[0])
        else:
            data.append(entry)

    return data

def abn_read_numpy(f, opts=None, first=0, last=-1):
    """Read a data block of numbers from a .abn file, and return it as a Numpy
    array.  The file pointer is placed after the last data element read."""

    try:
        import numpy as np
    except ImportError:
        print >> sys.stderr, "abn_read_numpy: could not load numpy module"
        return None

    # If f is a file object, fine.  Otherwise assume it is a filename, and open the file.
    if isinstance(f, file):
        filename = f.name
    else:
        filename = f
        f = open(filename, "rb")
    if not f:
        print >> sys.stderr, "abn_read: could not open '%s'" % filename
        return None

    # Read .abn header
    n = 0
    size = 0
    endian = ''
    fmt = ''
    line = f.readline()
    while line:
        m = re.match(begindata_pattern, line)
        if m:
            n = int(m.group(1))
            size = int(m.group(2))
            endian = m.group(3)
            fmt = m.group(4) 
            break
        elif opts is not None:
            opts.update(cfg_read_line(line))
        line = f.readline()

    if not line:        # reached end of file without finding a data header
        print >> sys.stderr, "abn_read_numpy: could not find data header in '%s'" % filename
        return None

    # Check the format: for now we only accept simple formats, since dealing
    # with Numpy record types is a little complicated
    try:
        if endian == 'L':
            dt = np.dtype('<' + fmt)
        else:
            dt = np.dtype('>' + fmt)
    except:
        print >> sys.stderr, "abn_read_numpy: invalid format '%d'" % fmt
        return None

    # Check for sensible limits on requested data
    if first < 0 or first > n-1:
        first = 0
    if last < 0  or last > n-1:
        last = n-1

    # Skip ahead to first requested entry
    f.seek(first*size, SEEK_CUR)

    # Read binary data
    data = np.fromfile(f, dtype=dt, count=last-first+1)
    return data
