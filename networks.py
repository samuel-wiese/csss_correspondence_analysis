import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import os
import csv

from typing import List


def create_region_species_network(countries: List[str], taxonomic_groups: List[str] = None, full_latin: bool = False)\
		-> None:
	"""
	Saves a region-species network from MOL data.

	Parameters
	----------
	countries : List[str]
		The list of countries considered.
	taxonomic_groups : List[str]
		The considered taxonomic groups; all by default.
	full_latin : bool
		Whether to take the full latin name for specifying physics.
	"""

	# Read MOL data
	edges = []
	all_regions = []
	all_species = []
	for country in countries:
		_, _, regions = os.walk("mol_data/" + country + "/").__next__()
		regions = ["#" + r[0:-4] for r in regions]
		all_regions.extend(regions)
		for region in regions:
			with open("mol_data/" + country + "/" + region[1:] + ".csv", "r") as f:
				reader = csv.reader(f, delimiter=",")
				for i, line in enumerate(reader):
					if i < 4:
						continue
					tax = line[3]
					if full_latin:
						spicy = line[0]
					else:
						spicy = line[0][0:line[0].index(" ")]
					if taxonomic_groups is None or tax in taxonomic_groups:
						edges.append([region, spicy])
						if spicy not in all_species:
							all_species.append(spicy)

	# Create and save the network
	g = nx.Graph()
	g.add_nodes_from(all_regions, type="region")
	g.add_nodes_from(all_species, type="spicy")
	g.add_edges_from(edges)

	nx.write_gpickle(g, "networks/region_species.pkl")


# TODO: need to fix this or restrict to specific types
def plot_region_species_network() -> None:
	"""
	Plots the region-species network.
	"""

	g = nx.read_gpickle("networks/region_species.pkl")
	regions, species = nx.bipartite.sets(g)
	pos = {}
	pos.update((node, (1, index / len(regions))) for index, node in enumerate(regions))
	pos.update((node, (2, index / len(species))) for index, node in enumerate(species))
	nx.draw(g, pos=pos)
	plt.show()


def plot_node_degrees() -> None:
	"""
	Plots node degree distribution of the region-species network.
	"""

	# Get data
	g = nx.read_gpickle("networks/region_species.pkl")
	species_degrees = [g.degree[x] for x, y in g.nodes(data=True) if y["type"] == "spicy"]

	# Save
	plt.figure(figsize=(7, 4))
	plt.hist(species_degrees, bins=50, color="k")
	plt.xlabel("Regions")
	plt.ylabel("Number of Species by the\nNumber of Regions that are in")
	plt.tight_layout()
	plt.savefig("img/node_degrees.pdf", dpi=300, bbox_inches="tight")
	plt.clf()


create_region_species_network(["United States of America"])
# plot_region_species_network()
plot_node_degrees()
