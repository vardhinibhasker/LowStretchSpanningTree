#Computing a Low Stretch Spanning Tree for an unweighted graph.
import networkx as nx
from math import log
import star

class Tree:
    """ Unweighted low stretch spanning tree."""

    def __init__(self, G, center, T):
        if center not in G:
            raise nx.NodeNotFound('The given center is not in G.')
        self.G = G
        self.center = center
        self.T = T

    def __str__(self):
        return f"Tree {self.T.edges}."

    @classmethod
    def create(cls, G, center):
        """ Creates our unweighted low stretch spanning tree from a given graph,
        G and a center node.
        """
        T = nx.Graph()
        alpha = 1/(2*log(len(G)+6, 4/3))
        edges = Tree.recursive_star(G, center, alpha)
        T.add_edges_from(edges)
        return cls(G, center, T)

    def recursive_star(G, center, alpha):
        """ Recursively applies a star decomposition onto a graph to split it up
        into star components (see Star class), then onto subgraphs formed by
        identified components until there are no cycles left. Edges between star
        components are known as bridges and form the edges of the final tree.

        Parameters:
            alpha : float. Multiplicative upper bound on the radius of the new
            tree against the radius of the original graph. [See Elkin et al.
            paper for description].

        Returns:
            bridges : tuple. Edges of the final tree.
        """
        num_edges = G.size()
        if num_edges <= 2:
            return list(G.edges)
        vertex_sets, bridges = star.star_decomp(G, center, 1/3, alpha, num_edges)
        for vs, c in vertex_sets:
            bridges += Tree.recursive_star(G.subgraph(vs), c, alpha)
        return bridges

    # TODO: currently in star
    def get_bridges(self):
        pass

    #TODO: Return star components to perform recursion on.
    def get_star_components(self):
        pass

    # TODO
    def visualise(self):
        pass

    # TODO
    def calculate_stretch(self):
        pass
