import heapq
import numpy as np
import msprime
import tskit

ts = msprime.simulate(10, recombination_rate=2, random_seed=1)
nodes = ts.tables.nodes
edges = ts.tables.edges
sample = [0,1,2]

ts1 = simplify(sample, nodes, edges, ts.sequence_length)

class Segment(object):
    def __init__(self, left, right, node):
        assert left < right
        self.left = left
        self.right = right
        self.node = node
    def __lt__(self, other): 
        return self.left < other.left
    def __repr__(self):
        return repr((self.left, self.right, self.node))

def simplify(S, Ni, Ei, L):
    # initialize tree sequence tables
    tables = tskit.TableCollection(L)
    No = tables.nodes
    Eo = tables.edges
    # initialize the segment list -> a list for each node in the current
    # (unsimplified) ts
    A = [[] for _ in range(len(Ni))]  
    Q = []
    for u in S:  # for each node in the extant sample
        # add the node to the node table
        v = No.add_row(time=Ni.time[u], flags=1)  
        # add the segment to the A list
        A[u] = [Segment(0, L, v)]
    for u in range(len(Ni)):  # for each node in the complete ts
        for e in [e for e in Ei if e.parent == u]:  # for each edge where the curr node is parent
            for x in A[e.child]:  # for each segment associated with the child
                if x.right > e.left and e.right > x.left:
                    # if x lies within the interval defined by e
                    y = Segment(max(x.left, e.left), min(x.right, e.right), x.node)
                    heapq.heappush(Q, y)
        v = -1
        while len(Q) > 0:  # exhaust the qeue
            l = Q[0].left
            r = L
            X = []
            while len(Q) > 0 and Q[0].left == l:
                x = heapq.heappop(Q)
                X.append(x)
                r = min(r, x.right)
            if len(Q) > 0:
                r = min(r, Q[0].left)
            if len(X) == 1:
                x = X[0]
                alpha = x
                if len(Q) > 0 and Q[0].left < x.right:
                    alpha = Segment(x.left, Q[0].left, x.node)
                    x.left = Q[0].left
                    heapq.heappush(Q, x)
            else:
                if v == -1:
                    v = No.add_row(time=Ni.time[u])
                alpha = Segment(l, r, v)
                for x in X:
                    Eo.add_row(l, r, v, x.node)
                    if x.right > r:
                        x.left = r
                        heapq.heappush(Q, x)
            A[u].append(alpha)
    E = list(Eo)
    Eo.clear()
    E.sort(key=lambda e: (e.parent, e.child, e.right, e.left))
    start = 0
    for j in range(1, len(E)):
        condition = (
                E[j-1].right != E[j].left or
                E[j-1].parent != E[j].parent or
                E[j-1].child != E[j].child)
        if condition:
            Eo.add_row(E[start].left, E[j-1].right, E[j-1].parent, E[j-1].child)
            start = j
    j = len(E)
    Eo.add_row(E[start].left, E[j-1].right, E[j-1].parent, E[j-1].child)
    return tables.tree_sequence()
