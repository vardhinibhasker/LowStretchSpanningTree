#Computing a Low Stretch Spanning Tree for an unweighted graph.
import networkx as nx
import matplotlib.pyplot as plt
import spalgs
from math import log, ceil, floor
import timeit
import numpy as np
import random
from tqdm import tqdm
from collections import defaultdict, deque

# G is the graph with the ball (i.e. original, no nodes removed)
# H is the graph where vertices will be removed
def get_ball(distances, radius):
    ball = []
    d = 0
    while d <= radius:
        ball+=distances[d]
        d+=1
    return ball

#grows ball -delta 1/3
def ball_cut(G, dist_from_center, rho, delta, edge_num, source):
    r = rho * delta
    c = log(edge_num+1, 2) / ((1- 2*delta)*rho)

    if r >= 1:
        ball = get_ball(dist_from_center, r)
        cut_size = nx.cut_size(G, ball)
        volume = len(G.edges(ball))
    else:
        ball = [source]
        cut_size = G.degree[source] #UNWEIGHTED
        volume = cut_size
    while cut_size > c*(volume+1):
        r += 1
        ball+=dist_from_center[floor(r)]
        cut_size = nx.cut_size(G, ball)
        volume = len(G.edges(ball))
    return r, ball

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


def get_cone(G, l, ideal):
    cone = ideal
    while l > 0:
        for item in list(nx.node_boundary(G, cone)):
            cone.append(item)
        l-=1
    return cone

def cone_props(G, cone, edge_num):
    cone_subgraph_size = G.subgraph(cone).size()
    cone_cut_size = nx.cut_size(G,cone)
    volume = cone_subgraph_size + cone_cut_size
    if cone_subgraph_size == 0: mu = (volume+1)*log(edge_num+1, 2)
    else: mu = (volume)*log(edge_num/cone_subgraph_size)
    return mu, cone_cut_size

def cone_cut(G, x, l, L, S, edge_num, distances):
    r = l
    ideal = get_ideal(G, S, x, distances)
    cone = ideal if r == 0 else get_cone(G, r, ideal)
    mu, cone_cut_size = cone_props(G, cone, edge_num)
    while cone_cut_size > mu/(L-l):
        for item in list(nx.node_boundary(G, cone)):
            cone.append(item)
        mu, cone_cut_size = cone_props(G, cone, edge_num)
        r+=1
    return r, cone

def S_neighbours(G, S):
    neighs = set()
    for node in S:
        neighs.update(G.adj[node])
    return neighs

def contracted_distances(G, S):
    neighs = S_neighbours(G, S)
    H = G.copy()
    H.remove_nodes_from(S)
    v = 's'
    H.add_node(v)
    for neigh in neighs:
        if neigh in S: continue
        else: H.add_edge(v, neigh)
    return nx.single_source_shortest_path_length(H, v)

def cone_decomp(H, S, Delta, edge_num):
    #G = H.copy()
    S_distances = contracted_distances(H, S)
    cones = []
    anchors = []
    while S:
        for node in S: anchor = node; break
        r, cone = cone_cut(H, anchor, 0, Delta, S, edge_num, S_distances)
        for node in cone:
            H.remove_node(node)
            if node in S: S.remove(node)
        cones += [(cone, anchor)]
        anchors.append(anchor)
    return cones, anchors

def get_bridges(G, center, anchors, cutoff):
    visited = {center}
    queue = deque([(center, 0, G.neighbors(center))])
    while queue:
        parent, depth_now, children = queue[0]     ## length of shortest path is depth_now
        try:
            child = next(children)
            if child in anchors:
                anchors.remove(child)
                yield (child, parent)
            if child not in visited:
                visited.add(child)
                if depth_now < cutoff:
                    queue.append((child, depth_now + 1, G.neighbors(child)))
        except StopIteration:
            queue.popleft()

def star_decomp(G, center, delta, eps, edge_num):
    H = G.copy()
    distances = spalgs.single_source_shortest_path_length(G, center)
    rad = max(distances)-1
    r0, ball = ball_cut(G, distances, rad, delta, edge_num, center)
    cutoff = floor(r0)
    S = set(nx.node_boundary(G, ball))
    H.remove_nodes_from(ball)
    cones, anchors = cone_decomp(H, S, eps*rad/2, edge_num)
    bridges = list(get_bridges(G, center, anchors, cutoff))
    partitions = [(ball, center)] + cones
    return partitions, bridges

def star_recurse(G, center, alpha):
    if G.number_of_nodes() <= 2:return list(G.edges)
    else:
        partitions, bridges = star_decomp(G, center, 1/3, alpha, G.size())
        for part in partitions:
            bridges += star_recurse(G.subgraph(part[0]), part[1], alpha)
        return bridges

#returns a list of edges
def unweighted_lst(G, center):
    node_num = len(G)
    alpha = 1/(2*log(node_num+6, 4/3))
    tree_edges = star_recurse(G, center, alpha)
    T = nx.Graph()
    T.add_nodes_from(G)
    T.add_edges_from(tree_edges)
    return T


def draw_star(G, partitions, node_num):
    #try boundarynorm instead?
    colour_list = ['r','g','b','c','m','y', 'w']
    node_colours = ['k' for i in range(node_num)]
    i=0
    print(partitions)
    for nodes, anc in partitions:
        colour = colour_list[i]
        i+=1
        for node in nodes:
            for node in nodes:
                node_colours[node] = colour
    nx.draw_networkx(G, node_color=node_colours)
    plt.show()


def stretch(T, u, v):
    return nx.shortest_path_length(T, source=u, target=v)

def ave_stretch(G, T):
    edge_num = G.size()
    sum = 0
    for u, v in G.edges:
        sum += nx.shortest_path_length(T, source=u, target=v)
    return (sum/edge_num)

def find_girth(G):
    seen = defaultdict(lambda: defaultdict(lambda: False))
    cc = nx.number_connected_components(G)
    girth = float('inf')
    for v in G.nodes:
        #if G.degree[v] == 1: continue #no cycle through v
        for u in G.adj[v]:
            if not seen[v][u]:
                #if G.degree[u] == 1: continue #no other path from v to u => no cycle
                G.remove_edge(u, v)
                try:
                    l = nx.shortest_path_length(G, v, u)
                except nx.NetworkXNoPath:
                    continue
                cycle_length = l + 1
                if cycle_length < girth: girth = cycle_length
                G.add_edge(u, v)
                seen[u][v] = True
                seen[v][u] = True
    return girth
