import networkx as nx
import random

def erdos_renyi(m, edge_prob):
    connected = False
    while not(connected):
        G = nx.erdos_renyi_graph(n, edge_prob)
        if nx.is_connected(G): connected=True
    return G

#  0 < edge_prob <= 1, ideally greater than 0.3 to increase likelhood of a connected graph.
# Weights uniformly distributed.
def weighted_erdos_renyi(max_weight, n, edge_prob):
    weights = [1/x for x in range(1, max_weight)]
    connected = False
    while not(connected):
        G = nx.erdos_renyi_graph(n, edge_prob)
        if nx.is_connected(G): connected=True
    for u,v in G.edges:
        G[u][v]["weight"] = random.choice(weights)
    return G
