import networkx as nx
import random

def weighted_random_regular(max_weight, degree, n):
    weights = [1/x for x in range(1, max_weight)]
    G = nx.random_regular_graph(degree, n)
    for u,v in G.edges:
        G[u][v]["weight"] = random.choice(weights)
    return G
