import networkx as nx

def small_test():
    G = nx.Graph()
    G.add_weighted_edges_from([(0, 1, 3.0), (1, 2, 7.5), (1, 3, 4.0)])
    return G
