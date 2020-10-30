#Computing a Low Stretch Spanning Tree for an unweighted graph.
import networkx as nx
from math import log
import matplotlib.pyplot as plt
import star

# returns tree edges
def star_recurse(G, center, alpha):
    num_edges = G.size()
    if num_edges <= 2: return list(G.edges)
    else:
        vertex_sets, bridges = star.star_decomp(G, center, 1/3, alpha, num_edges)
        for vs, c in vertex_sets:
            bridges += star_recurse(G.subgraph(vs), c, alpha)
        return bridges

#returns a low stretch spanning tree for unweighted graphs.
def unweighted_lst(G, center):
    if center not in G: raise nx.NodeNotFound('The center is not in G')
    num_nodes = len(G)
    alpha = 1/(2*log(num_nodes+6, 4/3))
    tree_edges = star_recurse(G, center, alpha)
    T = nx.Graph()
    T.add_nodes_from(G)
    T.add_edges_from(tree_edges)
    return T
