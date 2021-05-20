#Computing a Low Stretch Spanning Tree for an unweighted graph.
import networkx as nx
from math import log
import star

class Tree:
    """ Unweighted low stretch spanning tree."""
    def __init__(self, G, center):
        self.G = G
        self.center = center
        if center not in G:
            raise nx.NodeNotFound('The given center is not in G.')
        self.T = nx.Graph()

    def __str__(self):
        if self.T:
            return f"Tree {self.T.edges}."
        return f"Graph {self.G.edges}, with center {self.center}."

    def create_tree(self):
        """ Adds edges to our final unweighted low stretch spanning tree and
        returns it.

        Returns:
            T : nx.Graph.
        """
        alpha = 1/(2*log(len(self.G)+6, 4/3))
        edges = self.recursive_star(alpha)
        self.T.add_edges_from(edges)
        return self.T

    def recursive_star(self, alpha):
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
        num_edges = self.G.size()
        if num_edges <= 2:
            return list(self.G.edges)
        vertex_sets, bridges = star.star_decomp(self.G, self.center, 1/3, alpha, num_edges)
        for vs, c in vertex_sets:
            t = Tree(self.G.subgraph(vs), c)
            bridges += t.recursive_star(alpha)
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
