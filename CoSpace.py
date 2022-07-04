import networkx as nx
import numpy as np
import graph_tool.all as gt
from networks import create_region_species_network, take_monopartite_projections


class CoSpace:
    """

    Builds a co-occurence graph (or space)

    """

    def __init__(self, geo_groups, tax_groups):
        """

        :param geo_groups: The groups of regions considered (possible: Canada, China, Europe, United-States-of-America)
        :param tax_groups: The considered taxonomic groups; all if empty.
        """

        self.geo_groups = geo_groups
        self.tax_groups = tax_groups

        network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map, desc \
            = create_region_species_network(geo_groups, tax_groups)
        network_species, network_regions = take_monopartite_projections(network)
        self.network_species = network_species
        self.network_regions = network_regions

    def plot_species_space(self, edge_density_percentile: float):
        """

        Plots a space where vertices are species and edges are weighted by co-occurence between species in the same region

        :param edge_density_percentile: the percentile of edges to keep (ex: .1 keeps the top 10% of edges by weight)
        :return: Saves a pdf with an image of the network
        """
        name_save = f'./img/species_space_edp_{np.round(edge_density_percentile, decimals=2)}.pdf'
        self._plot_network(network=self.network_species, edge_density_percentile=edge_density_percentile, name_save=name_save, distance=10, vertex_size=0.2)

    def plot_region_space(self, edge_density_percentile: float):
        """
        Plots a space where vertices are regions and edges are weighted by co-occurence of species in the same region

        :param edge_density_percentile: the percentile of edges to keep (ex: .1 keeps the top 10% of edges by weight)

        :return: Saves a pdf with an image of the network
        """
        name_save = f'./img/region_space_edp_{np.round(edge_density_percentile, decimals=2)}.pdf'
        self._plot_network(network=self.network_regions, edge_density_percentile=edge_density_percentile, name_save=name_save, distance=8, vertex_size=3)

    @staticmethod
    def _nx_to_graph_tool(network: nx.Graph):
        """
        Transforms a networkx graph into a graph_tool graph
        :param network: a networkx graph
        :return: a graph_tool graph
        """

        # Build a graph tool graph
        adj = nx.to_numpy_matrix(network)
        idx = adj.nonzero()
        weights = adj[idx]
        g = gt.Graph(directed=False)
        g.add_edge_list(np.transpose(idx))
        gt.remove_parallel_edges(g)

        # add weights as an edge propetyMap
        ew = g.new_edge_property("double")
        ew.a = weights
        g.ep['edge_weight'] = ew

        # add node names
        names = np.array(list(network.nodes))
        node_names_vp = g.new_vertex_property("string", vals=names)
        g.vp['names'] = node_names_vp
        return g

    @staticmethod
    def _get_filtered_graph(network: nx.Graph, edge_density_percentile: float):
        """

        Filters out edges with small weight

        :param network: a weighted networkx graph
        :param edge_density_percentile: the percentile of edges to keep (ex: .1 keeps the top 10% of edges by weight)
        :return: a graph_tool graph whose edges are the top edge_density_percentile of edges by weight
        """

        g = CoSpace._nx_to_graph_tool(network=network)
        weights = np.array(g.ep['edge_weight'].a)

        edges_tree = np.array(gt.min_spanning_tree(g=g, weights=g.new_edge_property('float', vals=-1 * weights)).a)

        sorted_weights = np.sort(weights)[::-1]
        lower_bound = sorted_weights[int(edge_density_percentile * g.num_edges())]
        edges_strong = np.where(weights > lower_bound, True, False)

        edge_filter = np.logical_or(edges_tree, edges_strong)
        edge_filter_pm = g.new_edge_property('bool', vals=edge_filter)
        g.set_edge_filter(prop=edge_filter_pm)
        return g

    @staticmethod
    def _get_community_structure(g: gt.Graph):
        """

        Gets the community structre of the graph using the stochastic block model

        :param g: a graph
        :return: State with minimum description length for the stochastic block model
        """
        state = gt.minimize_blockmodel_dl(g)
        return state

    @staticmethod
    def _plot_network(network: nx.Graph, edge_density_percentile: float, name_save: str, distance: int, vertex_size: float):
        """

        Plots a network in two steps: (i) filters out edges with small weight (ii) finds community structure of the filtered network

        :param network: a graph
        :param edge_density_percentile: the percentile of edges to keep (ex: .1 keeps the top 10% of edges by weight)
        :param name_save: name with which you want to save the plot
        :param distance: the distance between the vertices
        :param vertex_size: the size of the vertices
        :return: Saves a pdf with an image of the network
        """
        g = CoSpace._get_filtered_graph(network=network, edge_density_percentile=edge_density_percentile)
        state = CoSpace._get_community_structure(g=g)
        pos = gt.fruchterman_reingold_layout(g, n_iter=10)
            # gt.sfdp_layout(g=g, eweight=g.ep['edge_weight'], groups=state.get_blocks(), K=distance)
        gt.graph_draw(g=g, pos=pos, output=name_save, vprops={'shape': state.get_blocks(), 'fill_color': state.get_blocks(), 'text': g.vp['names'], 'size': vertex_size})


if __name__ == '__main__':
    geo_groups = ["United-States-of-America", "China", "Canada", "Europe"]
    cp = CoSpace(geo_groups=geo_groups, tax_groups=[])
    # cp.plot_region_space(edge_density_percentile=.3)
    cp.plot_species_space(edge_density_percentile=.000000001)

