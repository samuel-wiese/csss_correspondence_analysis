import numpy as np
import pickle as pkl
import os
import csv

import networkx as nx
from networkx.algorithms import bipartite

from typing import List, Dict, Tuple


def create_region_species_network(groups: List[str], taxonomic_groups: List[str], full_latin: bool = False) \
		-> Tuple[nx.Graph, np.ndarray, np.ndarray, List[str], List[str], List[str], List[str], Dict, str]:
	"""
	Saves a region-species network from MOL data.

	Parameters
	----------
	groups : List[str]
		The groups of regions considered.
	taxonomic_groups : List[str]
		The considered taxonomic groups; all if empty.
	full_latin : bool
		Whether to take the full latin name for specifying physics.

	Returns
	-------
	network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, desc
		: Tuple[nx.Graph, np.ndarray, np.ndarray, List[str], List[str], List[str], List[str], Dict, str]
		The resulting region-species network in network and in (sorted) matrix form.
	"""

	# Check if we've done this already
	filename = groups[0]
	for region in groups[1:]:
		filename += "|" + region
	if len(taxonomic_groups) == 0:
		filename += "_all"
	else:
		filename += "_" + taxonomic_groups[0]
		for taxonomic_group in taxonomic_groups[1:]:
			filename += "|" + taxonomic_group
	network_filename = "networks/region_species_network_" + filename + ".pkl"
	matrix_filename = "matrices/region_species_network_" + filename + ".pkl"
	if os.path.exists(network_filename) and os.path.exists(matrix_filename):
		network = pkl.load(open(network_filename, "rb"))
		matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map\
			= pkl.load(open(matrix_filename, "rb"))
		return network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map,\
			   filename

	# Read MOL data
	edges = []
	all_regions = []
	all_species = []
	taxonomic_map = {}
	for group in groups:
		_, _, regions = os.walk("mol_data/" + group + "/").__next__()
		regions = ["#" + r[0:-4] for r in regions]

		# Removing Hawaii and Alaska from our analysis
		regions.remove("#Hawaii")
		regions.remove("#Alaska")

		all_regions.extend(regions)
		for region in regions:
			with open("mol_data/" + group + "/" + region[1:] + ".csv", "r") as f:
				reader = csv.reader(f, delimiter=",")
				for i, line in enumerate(reader):
					if i < 4:
						continue
					tax = line[3]
					if full_latin:
						spicy = line[0]
					else:
						spicy = line[0][0:line[0].index(" ")]
					if len(taxonomic_groups) == 0 or tax in taxonomic_groups:
						edges.append([region, spicy])
						if spicy not in all_species:
							all_species.append(spicy)
						if tax not in taxonomic_map:
							taxonomic_map[tax] = [spicy]
						else:
							taxonomic_map[tax].append(spicy)

	# Create the network and save it
	network = nx.Graph()
	network.add_nodes_from(all_regions, type="region")
	network.add_nodes_from(all_species, type="spicy")
	network.add_edges_from(edges)
	nx.write_gpickle(network, network_filename)

	# Computes the corresponding binary matrix
	regions = [x for x, y in network.nodes(data=True) if y["type"] == "region"]
	species = [x for x, y in network.nodes(data=True) if y["type"] == "spicy"]
	matrix = np.zeros((len(regions), len(species)))
	for ind_reg, region in enumerate(regions):
		for ind_spe, spicy in enumerate(species):
			if [region, spicy] in network.edges:
				matrix[ind_reg, ind_spe] = 1

	# Sort it
	matrix_sorted = np.zeros((len(regions), len(species)))
	regions_ind = np.argsort(np.sum(matrix, axis=1))[::-1]
	species_ind = np.argsort(np.sum(matrix, axis=0))[::-1]
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			matrix_sorted[i, j] = matrix[regions_ind[i], species_ind[j]]
	regions = np.array([r[1:] for r in regions])
	regions_sorted = regions[regions_ind]
	regions = list(regions)
	species_sorted = np.array(species)[species_ind]

	# Save it
	pkl.dump([matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map],
			 open(matrix_filename, "wb"))

	# Return it all
	return network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map, filename


def take_monopartite_projections(network: nx.Graph) -> Tuple[nx.Graph, nx.Graph]:
	"""
	Take monopartite projection of the complete bipartite region-species network.

	Parameters
	----------
	network : nx.Graph
		The region-species network.

	Returns
	-------
	network_species, network_regions : Tuple[nx.Graph, nx.Graph]
		Weighted monopartite projections.
	"""

	species = [x for x, y in network.nodes(data=True) if y["type"] == "spicy"]
	regions = [x for x, y in network.nodes(data=True) if y["type"] == "region"]
	network_species = bipartite.weighted_projected_graph(network, species, ratio=True)
	network_regions = bipartite.weighted_projected_graph(network, regions, ratio=True)

	return network_species, network_regions
