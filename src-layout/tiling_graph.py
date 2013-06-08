#!/usr/bin/env python

"""Plan:
* Read graph G
* Construct the grid assignment problem matrix
* Solve the assignment problem approximately to get a u->c assignment
* Create a routing graph: for every cell that is not filled replace it with a clique
* find a set of edge disjoint paths from {s_i} to {t_i} representing the unsatisfied edges

IDEAS:

y1. use sparse matrix somehow to speed up?
    - implictly store the matrix and do the multiplication
      on that implict matrix in the power method
    - separate matrix into (A + B + C) where A is sparse, and B 
      and C are easy to multiply
      
x2. weaken the requirement of adjacent => edge
    x look @ "loss" due to forbidden pairs
    - fan graph
4. limit # of routes passing through an edge
5. look for bug in hypercube(5)
6. local clean up / greedy improvement?
7. color by embedding into a line
8. offset routed edges (?)
9. handle weighted edges (?)
10. only iterate over "near" pairs of cells when doing distance matching

x3. add "border" before routing
x11. only add distance requirement for single hop edges
x12. color edges to avoid same color edges in a cell
x13. add constratints that put "close" nodes in close cells
"""

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import networkx as nx
import igraph
import itertools
import math
import sys
import colorsys
import random

#==========================================================================
# Represents a tiling object
#==========================================================================

class HexTiling(object):
    def __init__(self, x, y):
        """Create a tiling of at least x cells by y cells."""
        self.x = x
        self.y = y
        self.radius = 20
        self._init_hex_centers()

    def _init_hex_centers(self):
        """Initialize self.centers so that self.centers[c] = (x,y) center
        of cell c"""
        self.centers = {}
        c = x = y = 0
        for j in xrange(self.y):
            x = int(2 * self.hex_radius() + self.hex_side_len() / 2.0 if (j % 2 == 1) 
                        else self.hex_radius())
            y += int(self.hex_height() / 2);

            for i in xrange(self.x):
                self.centers[c] = (x,y)
                c += 1
                x += 3 * self.hex_radius()

    def hex_radius(self):
        """Radius of the circle enclosing the hexagon"""
        return self.radius

    def hex_height(self):
        """Height from side to side of the hexagon"""
        return 2 * self.hex_radius() * math.cos(2 * math.pi / 12)

    def hex_side_len(self):
        """The length of a tile side"""
        return 2 * self.hex_radius() * math.sin(2*math.pi / 12)

    def hex_sides(self):
        """Return the names of hex sides"""
        return ['se', 's', 'sw', 'nw', 'n', 'ne']

    def cell_center(self, c):
        """Return the center of cell c"""
        return self.centers[c]

    def cell_distance(self, c, d):
        """Return the distance between centers of c and d"""
        cx, cy = self.cell_center(c)
        dx, dy = self.cell_center(d)
        return math.sqrt((cx - dx)**2 + (cy - dy)**2)

    def neighbors(self, d):
        """Return dictionary of neighbors of d: D[side] = cell"""
        x = self.x
        if (d // x) % 2 == 1:
            D = {'nw': d - x,
                 'n' : d - 2*x,
                 'ne': d - x + 1,
                 'se': d + x + 1,
                 's' : d + 2*x,
                 'sw': d + x}
        else:
            D = {'nw': d - x - 1,
                 'n' : d - 2*x,
                 'ne': d - x,
                 'se': d + x,
                 's':  d + 2*x,
                 'sw': d + x - 1}

        # remove items that out out of range on top or bottom
        D = {x:p for x, p in D.iteritems() 
                if p >= 0 and p < self.num_cells()}

        # remove items that are out of range left or right
        first_in_this_row = x * ( d // x)
        if d == first_in_this_row and (d//x) % 2 == 0:
            if 'nw' in D: del D['nw']
            if 'sw' in D: del D['sw']
        last_in_this_row = first_in_this_row + self.x - 1
        if d == last_in_this_row and (d//x) % 2 == 1:
            if 'ne' in D: del D['ne']
            if 'se' in D: del D['se']
        return D

    def adjacent_cells(self):
        """Generator of pairs of adjacent cells"""
        return ( (c,d) for c in self.cells()
            for d in self.neighbors(c).itervalues()
                if c < d
        )

    def are_adjacent(self, c, d):
        """return true iff c and d share a side"""
        return d in set(self.neighbors(c).values())

    def cells(self):
        """Generator for all the cell ids"""
        return xrange(self.x * self.y)

    def num_cells(self):
        """The total number of cells"""
        return self.x * self.y

    def opposite_side(self, s):
        opp = {'n':'s', 's':'n', 'se':'nw', 'sw':'ne', 'ne':'sw', 'nw':'se'}
        return opp[s]

    def side_name(self, c, s):
        """Return the name of the side for cell c side s"""
        opp = {'n':'s', 's':'n', 'e':'w', 'w':'e'}
        N = self.neighbors(c)
        # if this is a side shared with a smaller cell N[s], use
        # that cell's name for this node
        if s in N and N[s] < c:
            c = N[s]
            s = "".join(opp[a] for a in s)
        return "s_%d_%s" % (c, s)

    def grid_coords(self, c):
        y = c // self.y
        x = c % self.y
        return (x,y)

    def coords_grid(self, x, y):
        assert 0 <= x < self.x
        assert 0 <= y < self.y
        c = y * self.x + x
        print "c=", c
        assert 0 <= c < self.num_cells()
        return c

#==========================================================================
# Assignment of nodes to cells via iterative linear assignment
#==========================================================================

def build_lap_graph(G, T, n2c, c2n, inertia):
    """Build a bipartite graph from nodes to cells"""
    
    n = G.number_of_nodes()
    LapGraph = igraph.Graph()
    LapGraph.add_vertices(n + T.num_cells())
    for u in G:
        LapGraph.vs[u]['type'] = True  # left half of bipartite graph
    for c in T.cells():
        LapGraph.vs[n+c]['type'] = False  # right half of bipartite graph

    for u in G:
        if random.random() < inertia:
            LapGraph.add_edge(u, n+n2c[u], weight=1000)
            continue

        for c in T.cells():
            cost_cu = 1
            bad_count = 0
            for d in T.neighbors(c).itervalues():
                if d in c2n and c2n[d] != u:
                    if G.has_edge(c2n[d], u):    
                        cost_cu += 1000
                    else:
                        bad_count += 1
                        #cost_cu = None
                        #break

            if cost_cu is not None and (bad_count == 0 or cost_cu > 1):
                for v in G.neighbors(u):
                    celldist = int(T.cell_distance(c, n2c[v]) / (2*T.hex_radius())) + 1
                    cost_cu += 500 / celldist

                assert cost_cu >= 0
                LapGraph.add_edge(u, n + c, weight=int(cost_cu))

    return LapGraph


def make_initial_assignments(G, T):
    """Construct some dummy initial assignments of n2c and c2n; not guaranteed to be good
    or even feasible."""
    # satisfy some edges
    n2c = {}
    c = 0
    for u,v in G.edges_iter():
        if u not in n2c and v not in n2c:
            n2c[u] = c
            n2c[v] = c+1
            c += 2

    # put unassinged nodes lined up at the end
    for u in G.nodes_iter():
        if u not in n2c:
            n2c[u] = c
            c += 1

    # construct c2n
    c2n = {c:u for u,c in n2c.iteritems()}

    assert len(set(u for u in n2c)) == len(G)

    return n2c, c2n


def find_blocked_sides(G, T, n2c, c2n):
    """find which sides of cells are adjacent to nodes with no edge"""
    Blocked = {}
    for u in G:
        c = n2c[u]
        Blocked[u] = set()
        for s, d in T.neighbors(c).iteritems():
            if d in c2n and not G.has_edge(u, c2n[d]):
                Blocked[u].add(s)
    return Blocked

                
def place_unplaced_nodes(G, T, n2c, c2n):
    """for every node in G that is not in n2c, put it someplace that is not forbidden"""

    # find the nodes that haven't been placed
    i = 0
    V = {}
    invV = []
    for u in G:
        if u not in n2c:
            V[u] = i
            invV.append(u)
            i += 1

    if len(V) == 0: return

    # build a bipartite graph
    L = igraph.Graph()
    L.add_vertices(len(V) + T.num_cells())

    for u in V:
        L.vs[V[u]]['type'] = 0
        for c in T.cells():
            # skip if c is already assigned
            if c in c2n: continue

            # skip if node is blocked
            for d in T.neighbors(c).itervalues():
                if d in c2n and not G.has_edge(c2n[d], u):
                    break
            else:
                L.add_edge(V[u], c + len(V), weight=1)
    for c in T.cells():
        L.vs[len(V) + c]['type'] = 1

    # find the matching
    M = L.maximum_bipartite_matching()
    for e in M.edges():
        u = invV[e.source]
        c = e.target - len(V)
        assert u not in n2c
        assert c not in c2n
        print "Placing: %d in cell %d" % (u,c)
        n2c[u] = c
        c2n[c] = u


def find_assignment_ilap(G, T, CUTOFF=0):
    """Use an iterative LAP apprach to place the nodes"""
    print >> sys.stderr, "Finding assignment via iterative LAP..."
    n = G.number_of_nodes()
    n2c, c2n = make_initial_assignments(G, T)
    prev_weight = 0 # some small number
    delta = 2*CUTOFF + 1
    print n, T.num_cells(), n*T.num_cells()
    inertia = 1
    global StatsIterations
    StatsIterations = 0
    while delta > CUTOFF:
        StatsIterations += 1
        
        L = build_lap_graph(G, T, n2c, c2n, inertia / 200.0)
        M = L.maximum_bipartite_matching(weights='weight')
        old_n2c = n2c.copy()
        n2c = {}
        c2n = {}
        w = 0.0
        for e in M.edges():
            u = e.source
            c = e.target - n
            assert c >= 0
            w += e["weight"]
            n2c[u] = c
            c2n[c] = u
        place_unplaced_nodes(G, T, n2c, c2n)
        assert len(n2c) == n
        assert len(c2n) == n
        changed = sum(1 for u in n2c if old_n2c[u] != n2c[u])
        print w, len(L.es), changed
        delta = abs(prev_weight - w)
        delta = changed
        prev_weight = w
        inertia += 0.5
    
    n2c, c2n = shift_layout(G, T, n2c)
    n2c, c2n = float_islands(G, T, n2c, c2n)
    Blocked = find_blocked_sides(G, T, n2c,c2n)
    return n2c, c2n, Blocked

def shift_layout(G, T, n2c):
    fail = False
    n2c_down = {}
    for u,c in n2c.iteritems():
        N = T.neighbors(c)
        if 's' not in N: 
            fail = True
            break
        N = T.neighbors(N['s'])
        if 'se' not in N: 
            fail = True
            break
        n2c_down[u] = N['se']

    # if we failed to shift it, we just dont' change anythign
    if fail:
        c2n = {c:u for u,c in n2c.iteritems()}
        return n2c, c2n
    else:
        c2n = {c:u for u,c in n2c_down.iteritems()}
        return n2c_down, c2n

#==========================================================================
# Greedy improvement by floating islands
#==========================================================================

def compute_islands(G, T, n2c, c2n):
    """Move each 'island' around to fit into a good place"""
    Islands = []
    Considered = set()
    for u,root in n2c.iteritems():
        if root in Considered: continue
        Considered.add(root)
        I = [root, (u,0,0)]
        rx,ry = T.grid_coords(root)
        queue = [root]
        while len(queue):
            c = queue[-1]
            queue = queue[:-1]
            N = T.neighbors(c)
            for s, d in N.iteritems():
                if d in Considered: continue
                if d in c2n and G.has_edge(c2n[c], c2n[d]):
                    dx,dy = T.grid_coords(d)
                    I.append( (c2n[d], dx-rx, dy-ry) )
                    queue.append(d)
                    Considered.add(d)
        Islands.append(I)

    print "Islands="
    for I in Islands:
        print '\t', I
    Nodes = set(u for I in Islands for u,x,y in I[1:]) 
    assert len(Nodes) == len(G)
    return Islands

        
def area(T, n2c):
    minx = maxx = miny = maxy = None
    for u,c in n2c.iteritems():
        x,y = T.grid_coords(c)
        minx = min(minx, x) if minx is not None else x
        maxx = max(maxx, x) if maxx is not None else x
        miny = min(miny, y) if miny is not None else x
        maxy = max(maxy, y) if maxy is not None else x

    # return the negative of the area
    return (maxx - minx) * (maxy - miny)

def edge_stretch(G, T, n2c, EndPoints=None):
    euk = 0.0
    for u,v in G.edges_iter():
        if EndPoints is None or u in EndPoints or v in EndPoints:
            xu,yu = T.grid_coords(n2c[u])
            xv,yv = T.grid_coords(n2c[v])
            euk += math.sqrt( (xu-xv)**2 + (yu-yv)**2 )
    return euk


def island_quality(I, cx,cy, G, T, n2c):
    # create a copy without the island placed
    n2c2 = dict(n2c) 
    for u,x,y in I[1:]: del n2c2[u]
    c2n = {d:u for u,d in n2c2.iteritems()}

    # place the island
    for u,x,y in I[1:]:
        # if the island would move ontop of another island
        if cx+x < 0 or cy+y < 0 or cx+x >= T.x or cy+y >= T.y: 
            return -2000000
        newcell = T.coords_grid(cx+x, cy+y)
        if newcell in c2n: return -30000000
        if newcell > T.num_cells(): return -40000000
        n2c2[u] = newcell

    return - area(T, n2c2)

    Island = set(u for u,c in I[1:])
    euk = edge_stretch(G, T, n2c2, Island)
    return -euk



def float_islands(G, T, n2c, c2n):
    Islands = compute_islands(G, T, n2c, c2n)
    moved = True
    i = 1
    while moved and i < 200:
        print "Island Pass", i
        moved = False
        for I in Islands:
            # find the best place for the island
            best_quality = None
            for x in xrange(T.x):
                for y in xrange(T.y):
                    q = island_quality(I, x,y, G, T, n2c) 
                    if q > best_quality or best_quality is None:
                        best_quality = q
                        best_loc = (x,y)
                print "Island::", I[1][0], I[0], best_loc, best_quality

            # move the island to the best place
            for u,x,y in I[1:]:
                n2c[u] = T.coords_grid(best_loc[0]+x,best_loc[1]+y)

            if best_loc != I[0]: 
                moved = True
                print "I=%d old=%s new=%s q=%f" % (I[1][0], str(I[0]), str(best_loc), best_quality)
            I[0] = best_loc
        i += 1

    c2n = {c:u for u,c in n2c.iteritems()}
    return n2c, c2n

#==========================================================================
# Assignment via spectral approximation to quadratic assignment
#==========================================================================

def print_assignmat(G, A):
    n,m = np.shape(A)
    assert n == m
    print "       ",
    for j in xrange(n):
        uj,cj = unpack_row_col(G, j)
        print "(%2d,%2d)" % (uj,cj),
    print
    for i in xrange(n):
        ui,ci = unpack_row_col(G, i)
        print "(%2d,%2d)" % (ui,ci),
        for j in xrange(n):
            print "%7g" % (A[i,j]),
        print
    print


def row_col(G, u, c):
    """assumes u and c are integers"""
    assert isinstance(c, (int, long))
    assert isinstance(u, (int, long))
    return G.number_of_nodes()*c + u


def unpack_row_col(G, i):
    # i = n*c + u
    c = i // G.number_of_nodes() 
    u = i % G.number_of_nodes()
    assert c * G.number_of_nodes() + u == i
    return u,c


def build_qap_matrix(G, T):
    """Return the QAP matrix as a sparse matrix. The idea is:
    if c and d are adjacent then u and v must be adjacent
        if u v not adj => c and d must not be adj
    if u and v are adjacent, we want c and d to be adjacent
    uv adj, cd adj: desired
    uv adj, cd not adj: not desired
    uv not adj, cd adj: forbidden
    uv not adj, cd not adj: fine
    """

    n = T.num_cells() * G.number_of_nodes()
    A = 2*np.ones((n,n))
    print "Graph n=", G.number_of_nodes()

    print >> sys.stderr, "Computing: Agreement between cells and graph"""
    SP_len = nx.all_pairs_dijkstra_path_length(G)
    #for u,v in itertools.product(G, G):
    #    if graphdist > 1: continue
    for u,v in G.edges_iter():
        graphdist = SP_len[u][v]
        for c,d in itertools.product(T.cells(), T.cells()):
            celldist = T.cell_distance(c, d) / T.hex_radius()
            q = abs(graphdist - celldist) + 1
            A[row_col(G,u,c),row_col(G,v,d)] = 2 + 100.0 / q
            A[row_col(G,v,d),row_col(G,u,c)] = 2 + 100.0 / q
            A[row_col(G,v,c),row_col(G,u,d)] = 2 + 100.0 / q
            A[row_col(G,u,d),row_col(G,v,c)] = 2 + 100.0 / q

    print >> sys.stderr, "Computing: forbidden assignments"
    for c,d in T.adjacent_cells():
        for u in G:
            for v in G:
                if u < v and not G.has_edge(u,v):
                    A[row_col(G,u,c),row_col(G,v,d)] = 0
                    A[row_col(G,v,d),row_col(G,u,c)] = 0
                    A[row_col(G,u,d),row_col(G,v,c)] = 0
                    A[row_col(G,v,c),row_col(G,u,d)] = 0

    print >> sys.stderr, "Computing: good assignments"
    for c,d in T.adjacent_cells():
        for u,v in G.edges_iter():
            A[row_col(G,u,c),row_col(G,v,d)] = 1000
            A[row_col(G,v,d),row_col(G,u,c)] = 1000
            A[row_col(G,u,d),row_col(G,v,c)] = 1000
            A[row_col(G,v,c),row_col(G,u,d)] = 1000
    assert (scipy.transpose(A) == A).all()
    assert (A >= 0.0).all()
    return A

def build_sparse_qap_matrix(G, T):
    """Return the QAP matrix as a sparse matrix"""
    n = T.num_cells() * G.number_of_nodes()
    for c,d in T.adjacent_cells():
        for u in G:
            for v in G:
                if u < v and not G.has_edge(u,v):
                    E.append( (row_col(G,u,c), row_col(G,v,d), 0) )
                    E.append( (row_col(G,v,d), row_col(G,u,c), 0) )
                    E.append( (row_col(G,u,d), row_col(G,v,c), 0) )
                    E.append( (row_col(G,v,c), row_col(G,u,d), 0) )
    
    for c,d in T.adjacent_cells():
        for u,v in G.edges_iter():
            E.append( (row_col(G,u,c), row_col(G,v,d), 1000) )
            E.append( (row_col(G,v,d), row_col(G,u,c), 1000) )
            E.append( (row_col(G,u,d), row_col(G,v,c), 1000) )
            E.append( (row_col(G,v,c), row_col(G,u,d), 1000) )

    # construct the sparse scipy matrix
    row = np.array([e[0] for e in E])
    col = np.array([e[1] for e in E])
    data = np.array([e[2] for e in E])
    #A = scipy.sparse.csr_matrix((data, (row,col)), shape=(n,n), dtype='d')


def leading_eig(A):
    """Get the sorted leading eigenvalue of A"""
    r,c = np.shape(A)
    assert r == c
    #w, v = scipy.sparse.linalg.eigsh(A, k=1, which='LA')
    #w, v = scipy.linalg.eig(A)
    w, v = scipy.linalg.eigh(A, eigvals=(r-1,r-1))
    mw = max(w)
    leading_eig = [j for j,x in enumerate(w) if x == mw][0]
    assert w[leading_eig] == max(w)
    if max(v[:,leading_eig]) < 0:
        v[:,leading_eig] *= -1
    assert min(v[:,leading_eig]) > 0

    return v[:, leading_eig]


def leading_eig_power(A, precision=1e-5):
    """Find the leading eigenvector using the power method"""

    # make an initial guess
    print >> sys.stderr, "Constructing initial guess"
    v = np.sum(A, 0)
    v /= np.linalg.norm(v)

    # temporary storage
    w = v.copy()

    # while we haven't converged
    print >> sys.stderr, "Using Power Method:",
    delta = 10000.0
    while delta > precision:
        # solve w = A*v and normalize w
        np.dot(A, v, out=w)
        nw = np.linalg.norm(w)
        w /= nw
        assert nw != 0

        # compute the change from last iteration
        delta = np.linalg.norm(v - w)
        v[:] = w
        print >> sys.stderr, delta,
    print >> sys.stderr
    return v


def round_assignment_relaxed(G, T, A, L):
    """Return an assignment of nodes to cells; allow some non-adjacent nodes to
    be assigned to adjacent cells."""

    LOSS_CUTOFF = 0.1

    FirstAssignment, Aunused, Bunused = round_assignment(G, T, A, L)

    Ldict = {unpack_row_col(G, i): x for i, x in L}

    AssignedNodes = {}
    AssignedCells = {}
    for i,x in L:
        u,c = unpack_row_col(G, i)

        # skip if one of these has already been matched
        if u in AssignedNodes or c in AssignedCells:
            continue

        # figure out which sides would block this assignment (if any)
        blockSet = set()
        for side, d in T.neighbors(c).iteritems():
            if d in AssignedCells and not G.has_edge(u, AssignedCells[d]):
                blockSet.add(side)

        # if there is something blocking, compute the loss
        loss = 0
        if len(blockSet) > 0:
            first_reward = Ldict[u, FirstAssignment[u]]
            assert first_reward >= 0 and x >= 0
            loss = (x - first_reward) / x

        # make this assignment if the cell is not blocked or the loss is too great
        print loss
        if len(blockSet) == 0 or loss >= LOSS_CUTOFF:
            AssignedNodes[u] = c
            AssignedCells[c] = u

    # post process to find blocked sides
    BlockedSides = {}
    for u,c in AssignedNodes.iteritems():
        BlockedSides[u] = set()
        for side, d in T.neighbors(c).iteritems():
            if d in AssignedCells and not G.has_edge(u, AssignedCells[d]):
                BlockedSides[u].add(side)
        
    assert len(AssignedNodes) == len(G)
    return AssignedNodes, AssignedCells, BlockedSides


def round_assignment(G, T, A, L):
    """Convert the fractional sorted list of assignments L to choices of
    assignments"""

    AssignedNodes = {}
    AssignedCells = {}
    BlockedSides = {}
    for i,x in L:
        u,c = unpack_row_col(G, i)

        # skip if one of these has already been matched
        if u in AssignedNodes or c in AssignedCells:
            continue

        # skip if there is a cell blocking this assignment
        for d in T.neighbors(c).itervalues():
            if d in AssignedCells and not G.has_edge(u, AssignedCells[d]):
                break
        else:
            AssignedNodes[u] = c
            AssignedCells[c] = u
            BlockedSides[u] = {}

    assert len(AssignedNodes) == len(G)
    return AssignedNodes, AssignedCells, BlockedSides


def find_assignment_qap(G, T):
    """Return an assignment of nodes to cells"""
    print >> sys.stderr, "Building QAP matrix..."
    A = build_qap_matrix(G, T)
    r,c = np.shape(A)
    print >> sys.stderr, "Solving quadratic assignment problem..."

    L = leading_eig_power(A)
    L = sorted(enumerate(L), key=lambda x: x[1].real, reverse=True)

    return round_assignment_relaxed(G, T, A, L)


#==========================================================================
# Routing of remaining edges
#==========================================================================

def construct_routing_graph(T, occupied):
    """
        occupied are cell ids that have a node assigned to them
        missing_edge are nodes that have some unhandled edge adj to them
    """
    NODE_CROSSING_COST = 10
    EDGE_COST = 1
    OPP = {'n':'s', 's':'n', 'e':'w', 'w':'e'}

    # construct the weighted clique graph
    H = nx.Graph()
    for c in T.cells():
        w = EDGE_COST + (NODE_CROSSING_COST if c in occupied else 0)
        for s1,s2 in itertools.combinations(T.hex_sides(),2):
            H.add_edge(T.side_name(c,s1), T.side_name(c,s2), weight=w)

        ##c for s,d in T.neighbors(c).iteritems():
        ##c    opp = "".join(OPP[a] for a in s)
        ##c    H.add_edge(T.side_name(c, s), T.side_name(d, opp), weight=0.0)

    return H


def construct_routing_graph_edges(T):
    H = nx.Graph()

    sides = T.hex_sides()
    assert len(sides) == 6

    for c in T.cells():
        # add the side-to-side edges
        for s1, s2 in itertools.combinations(sides, 2):
            H.add_edge("s_%d_%s" % (c,s1), "s_%d_%s" % (c,s2))
        
        for i in xrange(6):
            cs1 = sides[i % 6]
            cs2 = sides[(i + 1) % 6]
            corner = "c_%d_%s_%s" % (c,cs1, cs2)

            # add corner-to-side edges
            for s in sides:
                if s != cs1 and s != cs2:
                    H.add_edge(corner, "s_%d_%s" % (c, s))

            # get the cells across from this corner
            cn1 = T.neighbors(c)[cs1] if cs1 in T.neighbors(c) else None
            cn2 = T.neighbors(c)[cs2] if cs2 in T.neighbors(c) else None

            if cn1 is not None and cn2 is not None:
                junction = "j_%d_%d_%d" % tuple(sorted([c,cn1,cn2]))
            
                # add corner to junction
                H.add_edge(corner, junction)
            
            # add backward channel edges
            if cn1 is not None:
                cm, dm = min(c,cn1), max(c,cn1)
                g = "g_%d_%d" % (cm, dm)
                H.add_edge(junction, g)              # para channel
                H.add_edge("s_%d_%s" % (c, cs1), g)  # perp channel
            
            # add forward channel edges
            if cn2 is not None:
                cm, dm = min(c,cn2), max(c,cn2)
                g = "g_%d_%d" % (cm, dm)
                H.add_edge(junction, g)             # para channel
                H.add_edge("s_%d_%s" % (c, cs2), g) # para channel
            
    return H


def weight_routing_graph(RGraph, G, c2n):
    NODE_CROSSING_COST = 10
    EDGE_COST = 2
    SMALL_GAP_COST = 0.5
    LONG_GAP_COST = NODE_CROSSING_COST / 2.0

    ToDelete = []
    for u, v in RGraph.edges_iter():
        if u > v: u, v = v, u
        us = u.split('_')
        vs = v.split('_')
        etype = us[0] + vs[0]

        # side-to-side edges 
        if etype == 'ss':      
            cu, cv = int(us[1]), int(vs[1])
            assert cu == cv
            w = NODE_CROSSING_COST if cu in c2n else EDGE_COST

        # corner-to-side edges
        elif etype == 'cs':   
            cu, cv = int(us[1]), int(vs[1])
            assert cu == cv
            w = NODE_CROSSING_COST if cu in c2n else EDGE_COST

        # side-to-gap edges
        elif etype == 'gs':
            a, b = int(us[1]), int(us[2])
            c = int(vs[1])
            assert c == a or c == b
            # if cells a and c are blocked, add weight, else delete
            if a not in c2n or b not in c2n or G.has_edge(c2n[a], c2n[b]):
                w = None
            else:
                w = SMALL_GAP_COST

        elif etype == 'gj':   # gap-to-junction edges
            a, b = int(us[1]), int(us[2])
            if a not in c2n or b not in c2n or G.has_edge(c2n[a], c2n[b]):
                w = None
            else:
                w = LONG_GAP_COST

        # corner-to-junction edges
        elif etype == 'cj':   
            cu = int(us[1])
            for j in (int(x) for x in vs[1:]):
                if cu not in c2n or j not in c2n or G.has_edge(c2n[cu], c2n[j]):
                    w = 0.0
                    break
            else:
                w = SMALL_GAP_COST
            
        if w is not None:
            RGraph.edge[u][v]['weight'] = w
        else:
            ToDelete.append((u,v))

    for u,v in ToDelete:
        RGraph.remove_edge(u,v)
    return RGraph


def route_remaining_edges(G, T, n2c, c2n):
    """Return routes for each of the remaining edges"""
    #H = construct_routing_graph_edges(T) ##
    #weight_routing_graph(H, G, c2n) ##
    H = construct_routing_graph(T, c2n) ## clk: this one works
    nx.write_edgelist(H, "hex.graph")

    Routes = []
    for u, v in G.edges_iter():
        c = n2c[u]
        d = n2c[v]
        print c,d
        assert 0 <= c < T.num_cells()
        assert 0 <= d < T.num_cells()
        # find the combination of sides that gives the shortest path
        best = best_source = best_target = None
        for s1, s2 in itertools.product(T.hex_sides(),T.hex_sides()):
            source = T.side_name(c, s1)
            target = T.side_name(d, s2)
            if source in H and target in H:
                splen = nx.dijkstra_path_length(H, source, target)
                if splen < best or best is None:
                    best = splen
                    best_source = source
                    best_target = target
        P = nx.dijkstra_path(H, best_source, best_target)
        if len(P) < 2: 
            print "Bad Path!", P, c, d, source, target
            print best_source, best_target, best
            print best_source in H, best_target in H

        assert all(H.has_edge(P[i], P[i+1]) for i in xrange(len(P)-1))
        #for i in xrange(len(P)-1):
        #    H[ P[i] ][ P[i+1] ]['weight'] += 0.1
        Routes.append(P)
    return Routes

def route_remaining_edges_simple(G, T, n2c):
    """The original routing function --- not used now"""
    #for u,v in G.edges_iter():
    #    if T.are_adjacent(n2c[u], n2c[v]):
    #        print 'edge (%d,%d) at %d,%d good' % (u,v,n2c[u], n2c[v])


    if G.number_of_edges() == 0: return []

    H = construct_routing_graph(T, set(n2c.values()))
    SP = nx.all_pairs_dijkstra_path(H)
    SP_len = nx.all_pairs_dijkstra_path_length(H)
    nx.write_edgelist(H, "hex.graph")

    # for every remaining edge
    Routes = []
    for u,v in G.edges_iter():
        c = n2c[u]
        d = n2c[v]
        # find the combination of sides that gives the shortest path
        best = bestp = None
        for s1,s2 in itertools.product(T.hex_sides(),T.hex_sides()):
            source = T.side_name(c,s1)
            target = T.side_name(d,s2)

            if SP_len[source][target] < best or best is None:
                best = SP_len[source][target]
                bestp = SP[source][target]
        #print >>sys.stderr, "Route %d - %d (%g) %s" % (u, v, best, ",".join(bestp)) 
        Routes.append(bestp)
    return Routes

#==========================================================================
# Coloring of edges
#==========================================================================

def color_graph(G):
    """Return a color (integer) for every u in V to maximize so that
    the endpoints of no edge are the same color."""
    Nodes = sorted(G.nodes(), cmp=lambda a,b: cmp(G.degree(b), G.degree(a)))

    print >> sys.stderr, "%d degree=%d" % (Nodes[0], G.degree(Nodes[0]))
    
    k = 0
    palette = set()
    for u in Nodes:
        C = set(G.node[v]["color"] for v in G.neighbors(u) if "color" in G.node[v])
        # if there is no unused color
        Avail = palette - C
        if len(Avail) == 0:
            palette.add(k)
            G.node[u]["color"] = k
            k += 1
        else:
            G.node[u]["color"] = min(Avail)

    # return the number of colors used
    return k


def k_different_colors(k):
    """Return a list of k colors that are visually different from each other.  From:
    http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors"""

    colors = []
    for i in np.arange(0.0, 360.0, 360.0 / k):
        hue = i / 360.0
        lightness = (50 + np.random.rand() * 10) / 100.0
        saturation = (90 + np.random.rand() * 10) / 100.0
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append((255*rgb[0], 255*rgb[1], 255*rgb[2]))
        
    return colors
    
def color_routes(Routes):
    """Routes is a list of lists; each list is a path in the graph that will be
    drawn. Returns a list of the same length as Routes giving the colors."""
    if len(Routes) == 0: return []

    edge_name = lambda a,b: tuple(sorted((a, b)))

    ColorG = nx.Graph()

    # construct a dict mapping route to the edeges it contains
    Edges = {}
    for i, R in enumerate(Routes):
        Edges[i] = set(edge_name(R[j], R[j+1]) for j in xrange(len(R)-1))
        ColorG.add_node(i)

    Count = {}
    for E in Edges.itervalues():
        for e in E:
            if e not in Count: Count[e] = 0
            Count[e] += 1

    me = max(Count.itervalues())
    print >>sys.stderr, "Max degree edge=", me, ",".join(str(e) for e in Count if Count[e] == me)
    global StatsMaxEdgeDeg
    StatsMaxEdgeDeg = me

    for i, j in itertools.combinations(xrange(len(Routes)), 2):
        if len(Edges[i] & Edges[j]) > 0:
            ColorG.add_edge(i,j)

    nx.write_edgelist(ColorG, "routecolor.edg")

    k = color_graph(ColorG)
    print >> sys.stderr, "Required", k, "colors for routes"
    global StatColorsRequired
    StatColorsRequired = k
    
    Palette = k_different_colors(k)
    return [Palette[ColorG.node[i]["color"]] for i in xrange(len(Routes))]

#==========================================================================
# Computing of hex layout
#==========================================================================

def write_hex_layout(filename, G, T, node2cell, P, RouteColors, BlockedSides):
    """BlockedSides is a dictionary from node to a set of strings of the form nw, n, s, etc."""
    with open(filename, "wt") as out:
        print >> out, "T hex %d %d %d" % (T.x,T.y, G.number_of_nodes())

        # write the cell assignments
        for n,c in node2cell.iteritems():
            label = G.node[n]['label'] if 'label' in G.node[n] else str(n)
            label = label.replace(" ", "_")
            print >> out, "A %s %d %s %s" % (
                n, c, label, ",".join(BlockedSides[n])
                )

        # write the routes and their colors
        assert len(P) == len(RouteColors)
        for i, p in enumerate(P):
            if len(p) < 2: 
                print "WARNING: skipping an edge b/c there is no route"
                continue
            print >> out, "P c=%s %s" % (
                ",".join(str(x) for x in RouteColors[i]), " ".join(p)
            )


def hex_layout(G, x, y, filename):
    # construct a tiling
    n = G.number_of_nodes()
    m = G.number_of_edges()
    T = HexTiling(x,y)
    global StatsNumberEdges
    StatsNumberEdges = m

    # assign nodes to cells
    node2cell, cell2node, blockedSides = find_assignment_ilap(G, T)

    # route the rest of the edges

    print >>sys.stderr, "Started with %d edges." % (G.number_of_edges())
    # remove the edges from G that we have taken care of
    G.remove_edges_from([(u,v)
        for u,v in G.edges_iter()
            if T.are_adjacent(node2cell[u], node2cell[v])
        ])

    print >>sys.stderr, "%d edges remain." % (G.number_of_edges())

    global StatSatisfied
    StatSatisfied = m - G.number_of_edges()

    print >> sys.stderr, "Routing remaining edges..."
    P = route_remaining_edges(G, T, node2cell, cell2node)
    print "P=", P

    # choose colors for the edges
    RouteColors = color_routes(P)

    # write out the solution
    print >> sys.stderr, "Saving..."
    write_hex_layout(filename, G, T, node2cell, P, RouteColors, blockedSides)

#==========================================================================
# Main Program
#==========================================================================


def print_layout_statistics(G, T, n2c, c2n):
    print "Nodes:", G.number_of_nodes()
    print "Edges:", StatsNumberEdges
    print "Satisfied:", StatSatisfied
    print "Iterations:", StatsIterations
    print "ColorsRequired:", StatColorsRequired
    print "MaximumRouteDeg:", StatsMaxEdgeDeg
    print "Area:", area(T, n2c)
    #print "Crossings:",
    print "TotalEdgeLen:", edgeLen(G, T, n2c)

def number_nodes(G):
    H = nx.Graph()
    M = {u:i for i,u in enumerate(G)}
    for u in G:
        H.add_node(M[u], label=str(u))
    for u,v in G.edges_iter():
        H.add_edge(M[u], M[v])
    return H
        
def main():
    name = sys.argv[1]
    #graph=sys.argv[1]
    #G = nx.read_edgelist(graph)
    #G = nx.path_graph(10)
    #G = number_nodes(nx.hypercube_graph(4))
    #G = number_nodes(nx.grid_graph(dim=[4,4]))
    #G = number_nodes(nx.complete_graph(7))
    #G = number_nodes(nx.random_powerlaw_tree(300, tries=10000))
    #G = number_nodes(nx.read_gml(sys.argv[2]))
    #G = number_nodes(nx.florentine_families_graph())
    #G = number_nodes(nx.read_adjlist(sys.argv[2]))
    G = number_nodes(nx.read_edgelist(sys.argv[2]))
    nx.write_edgelist(G, name + "-in.graph")
    hex_layout(G, 15, 20, name + ".layout")

if __name__ == "__main__": main()
