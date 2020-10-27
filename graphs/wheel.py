import networkx as nx
import random

def weighted_wheel(max_weight, n):
    weights = [1/x for x in range(1, max_weight)]
    G=nx.wheel_graph(n)
    for u,v in G.edges:
        G[u][v]["weight"] = random.choice(weights)
    return G
