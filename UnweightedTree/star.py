import networkx as nx
from ball import ball_cut
from cone import cone_cut
from math import floor
from collections import deque

def distances_to_center(G, center):
    seen = set()
    level = 0
    dists = {}
    nextlevel = {center}

    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        vs = []
        for v in thislevel:
            if v not in seen:
                seen.add(v)
                nextlevel.update(G.adj[v])
                vs.append(v)
        if vs: dists[level] = vs
        level += 1
    return dists, level-2

def boundary_neighbors(G, node_boundary):
    neighbors = set()
    for node in node_boundary:
        neighbors.update(G.adj[node])
    return neighbors

# returns dictionary of shortest path lengths to the node boundary
def contracted_distances(G, node_boundary):
    neighs = boundary_neighbors(G, node_boundary)
    H = G.copy()
    H.remove_nodes_from(node_boundary)
    v = 's'
    H.add_node(v)
    for neigh in neighs:
        if neigh in node_boundary: continue
        else: H.add_edge(v, neigh)
    return nx.single_source_shortest_path_length(H, v)

# decomposes graph into cones
def cone_decomp(H, node_boundary, Delta, num_edges):
    node_boundary_distances = contracted_distances(H, node_boundary)
    cones, anchors = [],[]
    while node_boundary:
        for node in node_boundary: anchor = node; break
        r, cone = cone_cut(H, anchor, 0, Delta, node_boundary, num_edges, node_boundary_distances)
        for node in cone:
            H.remove_node(node)
            if node in node_boundary: node_boundary.remove(node)
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

def star_decomp(G, center, delta, eps, num_edges):
    H = G.copy()
    distances, radius = distances_to_center(G, center)
    ball_radius, ball = ball_cut(G, distances, radius, delta, num_edges, center)
    node_boundary = set(nx.node_boundary(G, ball))
    H.remove_nodes_from(ball)
    cones, anchors = cone_decomp(H, node_boundary, eps*radius/2, num_edges)
    bridges = list(get_bridges(G, center, anchors, floor(ball_radius)))
    partitions = [(ball, center)] + cones
    return partitions, bridges
