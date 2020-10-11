from collections import deque, defaultdict
import networkx as nx
from itertools import count
from heapq import heappush, heappop

def _single_shortest_path_length(adj, firstlevel, cutoff):
    seen = {}                  # level (number of hops) when seen in BFS
    level = 0                  # the current level
    nextlevel = firstlevel     # dict of nodes to check at next level

    while nextlevel and cutoff >= level:
        thislevel = nextlevel  # advance to next level
        nextlevel = {}         # and start a new list (fringe)
        vs = []
        for v in thislevel:
            if v not in seen:
                seen[v] = level  # set the level of vertex v
                nextlevel.update(adj[v])  # add neighbors of v
                vs.append(v)
        yield (level, vs)
        vs = []
        level += 1
    del seen


def single_source_shortest_path_length(G, source, cutoff=None):
    if source not in G:
        raise nx.NodeNotFound('Source {} is not in G'.format(source))
    if cutoff is None:
        cutoff = float('inf')
    nextlevel = {source: 1}
    return dict(_single_shortest_path_length(G.adj, nextlevel, cutoff))


def shortest_path_to_set(G, S, source):
    visited = {source}
    depth_limit = len(G)
    queue = deque([(source, 0, G.neighbors(source))])
    while queue:
        parent, depth_now, children = queue[0]     ## length of shortest path is depth_now
        if parent in S:
            return depth_now
        try:
            child = next(children)
            if child in S:
                return (depth_now+1)
            if child not in visited:
                visited.add(child)
                if depth_now < depth_limit:
                    queue.append((child, depth_now + 1, G.neighbors(child)))
        except StopIteration:
            queue.popleft()


def in_ideal(G, S, source, cutoff):
    """A fast BFS node generator"""
    G_adj = G.adj
    seen = set()
    nextlevel = {source}
    level = 0
    while nextlevel and cutoff > level:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v in S: return False
            if v not in seen:
                seen.add(v)
                nextlevel.update(G_adj[v])
        level+=1
    return True

def dijkstra(G, anchor):
    dist = {}  # dictionary of final distances
    seen = {}
    fringe = []
    seen[anchor] = 0
    heappush(fringe, (0, anchor))
    while fringe:
        (d, v) = heappop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        for u, e in G._adj[v].items():
            cost = G[u][v]["weight"]
            if cost is None:
                continue
            vu_dist = dist[v] + cost
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                #if u in closer to s, continue
                seen[u] = vu_dist
                heappush(fringe, (vu_dist, u))
    return dist

def to_cone(G, anc_dict, D, cutoff, ideal):
    dist = {}  # dictionary of final distances
    seen = {}
    fringe = []
    for elem in anc_dict:
        heappush(fringe, (anc_dict[elem], elem))
        seen[elem] = anc_dict[elem]
    while fringe:
        (d, v) = heappop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        yield v
        for u, e in G._adj[v].items():
            if u in ideal: continue
            cost = G[u][v]["weight"]
            if cost is None:
                continue
            try:
                D[v][u]
                vu_dist = dist[v]
            except KeyError:
                vu_dist = dist[v] + cost
            if vu_dist > cutoff: continue
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:','negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                heappush(fringe, (vu_dist, u))



def dijkstra(G, anchor):
    dist = {}  # dictionary of final distances
    seen = {}
    fringe = []
    seen[anchor] = 0
    heappush(fringe, (0, anchor))
    while fringe:
        (d, v) = heappop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        for u, e in G._adj[v].items():
            cost = G[u][v]["weight"]
            if cost is None:
                continue
            vu_dist = dist[v] + cost
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                #if u in closer to s, continue
                seen[u] = vu_dist
                heappush(fringe, (vu_dist, u))
    return dist



def in_ideal(G, node, S, anchor):
    G_adj = G._adj
    push = heappush
    pop = heappop

    dist = {}  # dictionary of final distances
    S_dist = float('inf')
    fringe = []
    push(fringe, (0, node))
    while fringe:
        (d, v) = pop(fringe)
        if d > S_dist: return False
        if v in dist: continue  # already searched this node.
        dist[v] = d
        if v == anchor: return True
        elif v in S: S_dist = d
        else:
            for u in G_adj[v]:
                cost = G[v][u]["weight"]
                if cost is None:
                    continue
                vu_dist = dist[v] + cost
                push(fringe, (vu_dist, u))
    return False



def dijkstra_to_set(G, source, S, pred=None, paths=None):
    G_succ = G._adj

    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    fringe = []
    if source not in G:
        raise nx.NodeNotFound("Source {} not in G".format(source))
    seen[source] = 0
    push(fringe, (0, 0, source))
    while fringe:
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if v in S:
            #break
            return d
        for u, e in G_succ[v].items():
            cost = G[v][u]["weight"]
            if cost is None:
                continue
            vu_dist = dist[v] + cost
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]
                if pred is not None:
                    pred[u] = [v]
            elif vu_dist == seen[u]:
                if pred is not None:
                    pred[u].append(v)
    return dist

def floyd_warshall_cycles(G, weight='weight'):
    dist = defaultdict(lambda: defaultdict(lambda: float('inf')))
    for u, v, d in G.edges(data=True):
        e_weight = d.get(weight, 1.0)
        dist[u][v] = min(e_weight, dist[u][v])
        dist[v][u] = min(e_weight, dist[v][u])
    for w in G:
        dist_w = dist[w]  # save recomputation
        for u in G:
            dist_u = dist[u]  # save recomputation
            for v in G:
                d = dist_u[w] + dist_w[v]
                if dist_u[v] > d:
                    dist_u[v] = d
                if u == v: print("The new weight is", dist_u[v])
    return dict(dist)
