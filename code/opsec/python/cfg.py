
import sys

def cfg_read_line(line):
    # Strip whitespace
    s = line.strip()

    # Ignore blank lines or lines that start with '#'
    if len(s) == 0 or s.startswith('#'):
        return {}

    # Split string on '=' character
    fields = s.split('=', 1)
    if len(fields) != 2:
        print >> sys.stderr, "cfg_read_line: '%s' is not a valid configuration line" % line
        return {}

    # Strip whitespace and/or quotes
    key = fields[0].strip()
    val = fields[1].strip()
    if len(val) >= 2 and ((val[0] == "'" and val[-1] == "'") or (val[0] == '"' and val[-1] == '"')):
        val = val[1:-1]

    return {key: val}


def cfg_read(f):
    """Reads configuration options from 'f', which may be either the name of a
    configuration file or an open file object.  Lines beginning with a '#'
    character are ignored.  Returns a dictionary."""

    if not isinstance(f, file):
        f = open(f, "r")
    if not f:
        print >> sys.stderr, "cfg_read: could not open '%s'" % f.name
        return {}

    opts = {}
    for line in f:
        opts.update(cfg_read_line(line))
    return opts
