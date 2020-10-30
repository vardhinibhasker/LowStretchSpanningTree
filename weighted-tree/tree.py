#Computing a Low Stretch Spanning Tree for a weighted graph.
import networkx as nx
import matplotlib.pyplot as plt
from math import log, ceil
import spalgs
from collections import defaultdict
import random
from heapq import heappush, heappop
from itertools import count
from tqdm import tqdm
import timeit
import quotient
import numpy as np

def graph_copy(G):
    H = nx.Graph()
    H.add_nodes_from(list(G.nodes))
    for u, v in G.edges:
        H.add_edge(u,v)
        H[u][v]["weight"] = G[u][v]["weight"]
    return H

def get_ball(distances, radius):
    ball = []
    ext = {}
    for key, val in distances.items():
        if val <= radius: ball.append(key)
        else: ext[key] = val
    return ball, ext

def ball_cut(G, distances, rho, delta, edge_num, source):
    r = rho * delta
    c = log(edge_num+1, 2) / ((1- 2*delta)*rho)

    ball, ext = get_ball(distances, r)
    cut_size = nx.cut_size(G, ball)
    volume = len(G.edges(ball))

    while cut_size > c*(volume+1):
        key = next(iter(ext))
        r = ext[key]
        ball.append(key)
        del ext[key]
        cut_size = nx.cut_size(G, ball)
        volume = len(G.edges(ball))
    return ball

def get_ideal(H, b_boundary, anchor):
    ideal = []
    dist_dict = nx.single_source_dijkstra_path_length(H, anchor) #spalgs.dijkstra(H, anchor)
    for node in dist_dict:
        dist_vS = spalgs.dijkstra_to_set(H, node, b_boundary)
        if dist_dict[node]==dist_vS:
            ideal.append(node)
    return ideal

def contracted_node_edges(G, S):
    s_edges = defaultdict(lambda : float('inf'))
    for node in S:
        for neigh in G.adj[node]:
            weight = s_edges[neigh]
            current_weight = G[node][neigh]["weight"]
            if current_weight < weight: s_edges[neigh] = current_weight
    return s_edges

def get_FS(H, distances):
    D = nx.DiGraph()
    for u, v in H.edges:
        dist_uS = distances[u]
        dist_vS = distances[v]
        weight = H[u][v]["weight"]
        if weight + dist_uS == dist_vS: D.add_edge(u,v)
        elif dist_vS + weight == dist_uS: D.add_edge(v, u)
    return D

def reachable_nodes(D, anchor):
    ideal = set()
    seen = set()
    nextlevel = {anchor}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                ideal.add(v)
                seen.add(v)
                nextlevel.update(D.neighbors(v))
    return ideal


def ideal_to_cone(G, r, ideal):
    cone = set(ideal)
    for node in ideal:
        length, path = nx.single_source_dijkstra(G, node, cutoff=r)
        cone.update(path)
    return list(cone)

#cone = ideal_to_cone(G, 3, ideal)

def set_dist(G, S):
    s_edges = contracted_node_edges(G, S)
    H = graph_copy(G)
    H.remove_nodes_from(S)
    v = 's'
    H.add_node(v)
    for node in s_edges:
        if node in S: continue
        else: H.add_edge(v, node, weight=s_edges[node])
    return H, v, nx.single_source_dijkstra_path_length(H, v)

def my_ideal(G, S, anchor, D, distances):
    D.add_node(anchor)
    for neigh in G.adj[anchor]:
        weight = G[anchor][neigh]["weight"]
        if neigh in distances:
            if distances[neigh] == weight: D.add_edge(anchor, neigh)
    return reachable_nodes(D, anchor)

def contract_ideal(H, ideal):
    contracted_ideal = defaultdict(lambda : float('inf'))
    for node in ideal:
        for neigh in H.adj[node]:
            if neigh in ideal: continue
            weight = contracted_ideal[neigh]
            current_weight = H[node][neigh]["weight"]
            if current_weight < weight: contracted_ideal[neigh] = current_weight
    return contracted_ideal

def get_cone(r, H, S, x, D, distances):
    ideal = my_ideal(H, S, x, D, distances)
    if r == 0: cone = ideal
    else:
        contr_ideal = contract_ideal(H, ideal)
        new_cones = spalgs.to_cone(H, contracted_ideal, D, r)
        cone = ideal.update(new_cones)
    return cone

def cone_props(G, cone, edge_num):
    vol_edges = G.subgraph(cone).size()
    vol_nodes = len(G.edges(cone))
    cut_size = 0
    for u, v in nx.edge_boundary(G, cone):
        cut_size+= 1/(G[u][v]["weight"])
    if vol_edges == 0: mu = (vol_nodes+1)*log(edge_num+1, 2)
    else: mu = (vol_nodes)*log(edge_num/vol_edges)
    return mu, cut_size

def grow_cone(H, cone):
    rad = float('inf')
    new_elems = set()
    for u, v in nx.edge_boundary(H, cone):
        dist = H[u][v]["weight"]
        if dist < rad:
            rad = dist
            new_elems = {u,v}
        if dist == rad: new_elems.update({u,v})
    cone.update(new_elems)
    return rad, cone

# # need full graph for this
def cone_cut(H, x, l, L, S, edge_num, D, distances):
    r = l
    cone = get_cone(r, H, S, x, D, distances)
    mu, cut_size = cone_props(H, cone, edge_num)
    while cut_size > mu/(L-l):
        rad, cone = grow_cone(H, cone)
        mu, cut_size = cone_props(H, cone, edge_num)
        r += rad
    return r, cone

# we remove from H, G is original
def cone_decomp(G, S, Delta, edge_num):
    H = graph_copy(G)
    C, v, distances = set_dist(G, S)
    if len(C)==1: return [([s], s) for s in S]
    D = get_FS(C, distances)
    cones = []
    while S:
        for node in S: anchor = node; break
        r, cone = cone_cut(H, anchor, 0, Delta, S, edge_num, D, distances)
        for node in cone:
            H.remove_node(node)
            D.remove_node(node)
            if node in S: S.remove(node)
        cones += [(cone, anchor)]
    return cones

def get_merged(neigh_dict):
    visited = set()
    partition = []
    for node in neigh_dict:
        if node not in visited:
            new_class = set()
            nodes = set([node])
            while nodes:
                node = nodes.pop()
                visited.add(node)
                nodes |= neigh_dict[node] - visited
                new_class.update([node])
            partition.append(frozenset(new_class))
    return partition

def get_contractions(G, const):
    neigh = defaultdict(set)
    for u, v in G.edges:
        if G[u][v]["weight"] < const:
            neigh[u].update([u,v]) ; neigh[v].update([u,v])
        else:
            neigh[u].update([u]) ; neigh[v].update([v])
    return get_merged(neigh)

def get_edge_weights(b, c):
    min_weight = float('inf')
    for u in b:
        for v in c:
            try:
                w = G[u][v]["weight"]
                if w < min_weight: min_weight=w
            except KeyError:
                continue
    if min_weight==float('inf'): print("BUG, edge weights not properly setting")
    return {'weight': min_weight}

def get_contracted_graph(G, const):
    contractions = get_contractions(G, const)
    Q = nx.quotient_graph(G, contractions, edge_data=get_edge_weights)
    renaming_dict = {}
    labels={}
    for count, contraction in enumerate(contractions):
        labels[contraction] = count
        for node in contraction:
            renaming_dict[node] = count
    Q = nx.relabel_nodes(Q, labels)
    return Q, renaming_dict, contractions

def quotient_g(G, renaming):
    edges = defaultdict(lambda: defaultdict(lambda: {'weight':float('inf')}))
    for u,v in G.edges:
        b = renaming[u]
        c = renaming[v]
        min_weight = edges[b][c]["weight"]
        w = G[u][v]["weight"]
        if w < min_weight: edges[b][c]["weight"] = w
    return nx.Graph(edges)

def contr_graph2(G, const):
    contractions = get_contractions(G, const)
    renaming_dict = {}
    for count, contraction in enumerate(contractions):
        for node in contraction:
            renaming_dict[node] = count
    Q=quotient_g(G, renaming_dict)
    return Q, renaming_dict, contractions

def update_partition(G, v_set, contractions):
    n_partition = []
    for v in v_set:
        n_partition += list(contractions[v])
    return n_partition

def quotient_star(Q, centerQ, delta, eps, rad, G, contractions):
    H = graph_copy(Q)
    pred, distances = nx.dijkstra_predecessor_and_distance(Q, centerQ)
    edge_num = Q.size()
    ball = ball_cut(Q, distances, rad, delta, edge_num, centerQ)
    S = set(nx.node_boundary(Q, ball))
    H.remove_nodes_from(ball)
    cones = cone_decomp(H, S, eps*rad/2, edge_num)
    return pred, cones, ball

def update_bridge(G, anc, b, contractions):
    new_edges = []
    ball = contractions[b]
    cone = contractions[anc]
    for u, v, d in G.edges(ball | cone, data=True):
        if u in ball and v in cone:
            new_edges.append((u, v, d.get('weight', 1)))
        if v in ball and u in cone:
            new_edges.append((v, u, d.get('weight', 1)))
    x, y, z = min(new_edges, key = lambda t : t[2])
    return (x, y, {"weight" : z})

def new_br2(G, anc, b, contractions, weight):
    ball = contractions[b]
    cone = contractions[anc]
    for u, v, d in G.edges(ball | cone, data=True):
        if u in ball and v in cone:
            if d.get('weight',1) == weight: return (u, v, {"weight" : weight})
        if v in ball and u in cone:
            if d.get('weight',1) == weight: return (v, u, {"weight" : weight})

def star_decomp_contracted(G, Q, centerQ, delta, eps, contractions, rad, renaming):
    predQ, conesQ, ballQ = quotient_star(Q, centerQ, delta, eps, rad, G, contractions)
    bridges = []
    partitions = []
    for v_set, anc in conesQ:
        v = predQ[anc][0]
        bridge = new_br2(G, anc, v, contractions, Q[anc][v]["weight"])
        bridges.append(bridge)
        new_cone = update_partition(G, v_set, contractions)
        partitions+=[(new_cone, bridge[1])]
    ball = update_partition(G, ballQ, contractions)
    return ball, partitions, bridges

def G_radius(G, center):
    radius = 0
    dist = {}  # dictionary of final distances
    seen = {}
    fringe = []
    seen[center] = 0
    heappush(fringe, (0, center))
    while fringe:
        (d, v) = heappop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if d > radius: radius=d
        for u, e in G._adj[v].items():
            vu_weight = e.get('weight', 1)
            if vu_weight is None:
                continue
            vu_dist = d + vu_weight
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                heappush(fringe, (vu_dist, u))
    return radius

def star_recurse(G, center, beta):
    node_num = len(G)
    if node_num <= 2: return list(G.edges(data=True))
    else:
        rad = G_radius(G, center)
        const = (beta*rad)/node_num
        Q, renaming, contractions = contr_graph2(G, const)
        new_center = renaming[center]
        ball, partitions, bridges = star_decomp_contracted(G, Q, new_center, 1/3, beta, contractions, rad, renaming)
        partitions += [(ball, center)]
        T = bridges
        for part, anc in partitions:
            subgr = G.subgraph(part)
            T += star_recurse(subgr, anc, beta)
        return T

#returns a list of edges
def weighted_lst(G, center):
    T = nx.Graph()
    node_num = len(G)
    beta = 1/(2*log(node_num+32, 4/3))
    tree_edges = star_recurse(G, center, beta)
    T.add_nodes_from(G)
    T.add_edges_from(tree_edges)
    for u, v in T.edges:
        T[u][v]["weight"] = G[u][v]["weight"]
    return T
