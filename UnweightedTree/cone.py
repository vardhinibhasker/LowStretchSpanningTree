import networkx as nx
from math import log

# outer base component of star decomposition
def get_ideal(G, S, anchor, distances):
    ideal = [anchor]
    # BFS from anchor
    seen = {anchor : 0}                  # level (number of hops) when seen in BFS
    level = 0                  # the current level
    nextlevel = set(G.adj[anchor])     # dict of nodes to check at next level
    temp=set()

    while nextlevel:
        thislevel = nextlevel  # advance to next level
        nextlevel = set()         # and start a new list (fringe)
        for v in thislevel:
            if v not in seen:
                if v in S: continue
                seen[v] = level  # set the level of vertex v
                d = distances[v]
                if d >= level:
                    ideal.append(v)
                    nextlevel.update(G.adj[v])  # add neighbors of v
                else: temp.update(G.adj[v])
            for u in temp:
                if u not in seen:
                    seen[u] = level+1
            temp = set()
        level += 1
    return ideal

# grows the ideal by radius l
def get_cone(G, l, ideal):
    cone = ideal
    while l > 0:
        for item in list(nx.node_boundary(G, cone)):
            cone.append(item)
        l-=1
    return cone

def cone_properties(G, cone, num_edges):
    cone_subgraph_size = G.subgraph(cone).size()
    cone_cut_size = nx.cut_size(G,cone)
    volume = cone_subgraph_size + cone_cut_size

    if cone_subgraph_size == 0: mu = (volume+1)*log(num_edges+1, 2)
    else: mu = (volume)*log(num_edges/cone_subgraph_size)
    return mu, cone_cut_size

#final cone
def cone_cut(G, x, l, L, S, num_edges, distances):
    r = l
    ideal = get_ideal(G, S, x, distances)
    cone = ideal if r == 0 else get_cone(G, r, ideal)
    mu, cone_cut_size = cone_properties(G, cone, num_edges)
    while cone_cut_size > mu/(L-l):
        for item in list(nx.node_boundary(G, cone)):
            cone.append(item)
        mu, cone_cut_size = cone_properties(G, cone, num_edges)
        r+=1
    return r, cone
