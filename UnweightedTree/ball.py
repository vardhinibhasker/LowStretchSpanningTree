import networkx as nx
from math import log

# central component of star decomposition
def get_ball(distances, radius):
    ball = []
    d = 0
    while d <= radius:
        ball+=distances[d]
        d+=1
    return ball

#grows ball
def ball_cut(G, dist_from_center, rho, delta, num_edges, source):
    radius = rho * delta
    c = log(num_edges+1, 2) / ((1- 2*delta)*rho)

    if radius >= 1:
        ball = get_ball(dist_from_center, radius)
        cut_size = nx.cut_size(G, ball)
        volume = len(G.edges(ball))
    else:
        ball = [source]
        cut_size = G.degree[source] #UNWEIGHTED
        volume = cut_size

    while cut_size > c*(volume+1):
        radius += 1
        ball+=dist_from_center[floor(radius)]
        cut_size = nx.cut_size(G, ball)
        volume = len(G.edges(ball))
    return radius, ball
