#!/usr/bin/env python

import re

def read_data ( fname ):
    seed = int(re.search ( r"data(\d+)", fname ).group(1))
    data = [ [float(x.split()[0]),int(x.split()[1]),int(x.split()[2])] for x in open(fname).readlines() ]
    return data,seed

if __name__ == "__main__":
    import sys
    d,s = read_data ( sys.argv[1] )
    print "seed:", s
    print "number of stimulus intensities:", len(d)
