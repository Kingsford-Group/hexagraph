#!/usr/bin/env python

import sys

var = sys.argv[1]
L = []
for f in sys.argv[2:]:
    with open(f) as inp:
        for line in inp:
            if line.startswith('S ' + var + ':'):
                print line.strip()
                s = line.strip().split()
                L.append(float(s[2]))

print var, sum(L) / float(len(L))

