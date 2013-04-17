#!/usr/bin/env python

"""Plan:

- color by embedding into a line

x Read graph G
x Construct the grid assignment problem matrix:
    minimize 
    A[(u,c),(v,d)] = INF for all cases where uv are not adjacent but cd are
    A[(u,c),(v,d)] = -1 if u and v are adjacent, and c and d are adjacent
* Solve the assignment problem approximately to get a u->c assignment

* Create a routing graph: for every cell that is not filled replace it with a clique

* find a set of edge disjoint paths from {s_i} to {t_i} representing the unsatisfied edges

if c and d are adjacent then u and v must be adjacent
    if u v not adj => c and d must not be adj

if u and v are adjacent, we want c and d to be adjacent

uv adj, cd adj: desired
uv adj, cd not adj: not desired
uv not adj, cd adj: forbidden
uv not adj, cd not adj: fine

IDEAS:

1. use sparse matrix somehow to speed up?
2. weaken the requirement of adjacent => edge
    - look @ "loss" due to forbidden pairs
    - fan graph
4. look for bug in hypercube(5)
5. local clean up / greedy improvement?
6. only add distance requirement for single hop edges

7. offset routed edges (?)
8. handle weighted edges (?)

x3. color edges to avoid same color edges in a cell
x5. add constratints that put "close" nodes in close cells

"""

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import networkx as nx
import itertools
import math
import sys
import colorsys

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
        return set(('nw', 'n', 'ne', 'se', 's', 'sw'))

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

    def side_name(self, c, s):
        """Return the name of the side for cell c side s"""
        opp = {'n':'s', 's':'n', 'e':'w', 'w':'e'}
        N = self.neighbors(c)
        # if this is a side shared with a smaller cell N[s], use
        # that cell's name for this node
        if s in N and N[s] < c:
            c = N[s]
            s = "".join(opp[a] for a in s)
        return "%s_%s" % (str(c), s)

    def grid_coords(self, c):
        y = c // self.y
        x = c % self.y
        return (x,y)

#==========================================================================
# Assignment of nodes to cells
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
    """Return the QAP matrix as a sparse matrix"""

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
            A[row_col(G,u,c),row_col(G,v,d)] = 2 + 10.0 / q
            A[row_col(G,v,d),row_col(G,u,c)] = 2 + 10.0 / q
            A[row_col(G,v,c),row_col(G,u,d)] = 2 + 10.0 / q
            A[row_col(G,u,d),row_col(G,v,c)] = 2 + 10.0 / q

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


def solve_qap(G, T, A):
    """Return an assignment of nodes to cells"""

    L = leading_eig(A)
    L = sorted(enumerate(L), key=lambda x: x[1].real, reverse=True)

    AssignedNodes = {}
    AssignedCells = {}
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
    assert len(AssignedNodes) == len(G)
    return AssignedNodes


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

    # construct the weighted clique graph
    H = nx.Graph()
    for c in T.cells():
        w = EDGE_COST + (NODE_CROSSING_COST if c in occupied else 0)
        for s1,s2 in itertools.combinations(T.hex_sides(),2):
            H.add_edge(T.side_name(c,s1), T.side_name(c,s2), weight=w)

    return H


def route_remaining_edges(G, T, n2c):
    #for u,v in G.edges_iter():
    #    if T.are_adjacent(n2c[u], n2c[v]):
    #        print 'edge (%d,%d) at %d,%d good' % (u,v,n2c[u], n2c[v])

    print >>sys.stderr, "Started with %d edges." % (G.number_of_edges())

    # remove the edges from G that we have taken care of
    G.remove_edges_from([(u,v)
        for u,v in G.edges_iter()
            if T.are_adjacent(n2c[u], n2c[v])
        ])

    print >>sys.stderr, "%d edges remain." % (G.number_of_edges())

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

    for i, j in itertools.combinations(xrange(len(Routes)), 2):
        if len(Edges[i] & Edges[j]) > 0:
            ColorG.add_edge(i,j)

    nx.write_edgelist(ColorG, "routecolor.edg")

    k = color_graph(ColorG)
    print >> sys.stderr, "Required", k, "colors for routes"
    
    Palette = k_different_colors(k)
    return [Palette[ColorG.node[i]["color"]] for i in xrange(len(Routes))]

#==========================================================================
# Computing of hex layout
#==========================================================================

def write_hex_layout(filename, G, T, node2cell, P, RouteColors):

    with open(filename, "wt") as out:
        print >> out, "T hex %d %d %d" % (T.x,T.y, G.number_of_nodes())

        # write the cell assignments
        for n,c in node2cell.iteritems():
            print >> out, "A %s %d %s" % (
                n, c, G.node[n]['label'] if 'label' in G.node[n] else str(n)
                )

        # write the routes and their colors
        assert len(P) == len(RouteColors)
        for i, p in enumerate(P):
            print >> out, "P c=%s %s" % (",".join(str(x) for x in RouteColors[i]), " ".join(p))


def hex_layout(G, x, y, filename):
    # construct a tiling
    n = G.number_of_nodes()
    T = HexTiling(x,y)

    # assign nodes to cells
    print >> sys.stderr, "Building QAP matrix..."
    A = build_qap_matrix(G, T)
    r,c = np.shape(A)
    print >> sys.stderr, "Solving quadratic assignment problem..."
    node2cell = solve_qap(G, T, A)

    # route the rest of the edges
    print >> sys.stderr, "Routing remaining edges..."
    P = route_remaining_edges(G, T, node2cell)

    # choose colors for the edges
    RouteColors = color_routes(P)

    # write out the solution
    print >> sys.stderr, "Saving..."
    write_hex_layout(filename, G, T, node2cell, P, RouteColors)

#==========================================================================
# Main Program
#==========================================================================

def number_nodes(G):
    H = nx.Graph()
    M = {u:i for i,u in enumerate(G)}
    for u in G:
        H.add_node(M[u], label=str(u))
    for u,v in G.edges_iter():
        H.add_edge(M[u], M[v])
    return H
        
def main():
    #graph=sys.argv[1]
    #G = nx.read_edgelist(graph)
    #G = nx.path_graph(10)
    #G = number_nodes(nx.hypercube_graph(4))
    #G = number_nodes(nx.grid_graph(dim=[4,4]))
    #G = number_nodes(nx.ladder_graph(20))
    G = number_nodes(nx.read_gml(sys.argv[1]))
    hex_layout(G, 12, 15, "hex.layout")

if __name__ == "__main__": main()
