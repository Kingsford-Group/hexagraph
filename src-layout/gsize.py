#!/usr/bin/env python

import sys
import networkx

G = networkx.read_edgelist(sys.argv[1])
print G.number_of_nodes(), G.number_of_edges()

