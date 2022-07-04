import numpy as np
import pandas as pd
import networkx as nx
import plotly.express as px
import csv

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from typing import List, Dict


def draw_binary_matrix(matrix: np.ndarray, desc: str) -> None:
	"""
	Draws the (sorted) matrix of species living in which regions.

	Parameters
	----------
	matrix : np.ndarray
		The (sorted) matrix.
	desc : str
		The description of what we're looking at.
	"""

	plt.xlabel("Species sorted by in how many regions they are")
	plt.ylabel("Regions sorted by how\nmany species they contain")
	plt.imshow(matrix, cmap=ListedColormap(['w', 'k'], N=2), vmin=0, vmax=1, aspect=0.5 * matrix.shape[1]
																					/ matrix.shape[0])
	plt.tight_layout()
	plt.savefig("img/sorted_matrix_" + desc + ".pdf", dpi=300, bbox_inches="tight")
	plt.clf()


def plot_node_degrees(network: nx.Graph, desc: str) -> None:
	"""
	Plots node degree distribution of the region-species network for both types of nodes.

	Parameters
	----------
	network: nx.Graph
		The network.
	desc : str
		The description of what we're looking at.
	"""

	# For regions
	region_degrees = [network.degree[x] for x, y in network.nodes(data=True) if y["type"] == "region"]
	plt.figure(figsize=(7, 4))
	plt.hist(region_degrees, bins=10, color="k")
	plt.xlabel("Species")
	plt.ylabel("Number of Regions by the\nNumber of Species they contain")
	plt.tight_layout()
	plt.savefig("img/region_degrees_" + desc + ".pdf", dpi=300, bbox_inches="tight")
	plt.clf()

	# For species
	species_degrees = [network.degree[x] for x, y in network.nodes(data=True) if y["type"] == "spicy"]
	plt.figure(figsize=(7, 4))
	plt.hist(species_degrees, bins=10, color="k")
	plt.xlabel("Regions")
	plt.ylabel("Number of Species by the\nNumber of Regions they live in")
	plt.tight_layout()
	plt.savefig("img/species_degrees_" + desc + ".pdf", dpi=300, bbox_inches="tight")
	plt.clf()


def print_metric_by_species(cou_s: np.ndarray, eci_s: np.ndarray, fit_s: np.ndarray, species: List[str],
							taxonomic_map: Dict) -> None:
	"""
	Print metric info on species pickyness.

	Parameters
	----------
	cou_s: np.ndarray
		The total number of regions the species doesn't live in.
	eci_s: np.ndarray
		The ECI-pickyness of each species.
	fit_s: np.ndarray
		The FC-pickyness of each species.
	species: List[str]
		The list of species.
	taxonomic_map : Dict
		The species in each taxonomic group.
	"""

	# Print some info about the metrics
	metrics = [cou_s, eci_s, fit_s]
	desc = ["Number of non-inhabitant Regions", "ECI", "FC"]
	species = np.array(species)
	print("")
	for i, metric in enumerate(metrics):
		print("For Pickyness - " + desc[i] + ":")

		# The most and least picky species
		print("- most picky species: " + ", ".join(species[np.argsort(metric)[::-1][0:5]]))
		print("- least picky species: " + ", ".join(species[np.argsort(metric)[0:5]]))

		# The most and least picky taxonomic groups, we're simply taking the mean since we have no weights sadly
		mean_pickyness = []
		for tax_group in taxonomic_map.keys():
			mean_pickyness.append(sum(metric[np.where(species == spicy)[0]] for spicy in taxonomic_map[tax_group])[0]
								  / len(taxonomic_map[tax_group]))
		print("- most picky taxonomic groups: "
			  + ", ".join(np.array(list(taxonomic_map.keys()))[np.argsort(mean_pickyness)[::-1][0:5]]))
		print("- least picky taxonomic groups: "
			  + ", ".join(np.array(list(taxonomic_map.keys()))[np.argsort(mean_pickyness)[0:5]]))

		print("\n")


def plot_metric_by_us_state(metric: np.ndarray, regions: List[str], desc: str, metric_desc: str,
							metric_desc_short: str) -> None:
	"""
	Plots some metric by US state.

	Parameters
	----------
	metric : np.ndarray
		List of the metric values.
	regions : List[str]
		List of region names.
	desc : str
		Description of current settings.
	metric_desc : str
		Description of the metric.
	metric_desc_short : str
		Short description of the metric.
	"""

	assert metric_desc is not None

	# Need to get US State abbreviations
	map_to_abbreviation = {}
	with open('trivia/us-states-territories.csv') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=',')
		for row in csv_reader:
			if row[0] != "Type":
				map_to_abbreviation[row[1]] = row[2]
	abb_regions = [map_to_abbreviation[r] for r in regions]

	# Draw it
	df = pd.DataFrame(list(zip(abb_regions, metric)), columns=['Region', "Metric"])
	fig = px.choropleth(df, locations='Region', locationmode="USA-states", scope="usa", color="Metric",
						color_continuous_scale="Viridis_r", hover_name="Region")
	fig.update_layout(margin=dict(l=0, r=0, b=0, t=0), autosize=True)
	fig.update_coloraxes(showscale=False)
	fig.write_image("img/map_usa_" + metric_desc_short + "_" + desc[desc.index("_") + 1:] + ".pdf")


# TODO: someone needs to finish this
def plot_clusters(network: nx.Graph, label: str, desc: str) -> None:
	"""
	Plots clusters for species and for regions to understand their similarity.

	Parameters
	----------
	network : nx.Graph
		A monopartite graph with edge weights.
	label : str
		Whether we're plotting species or regions.
	desc : str
		Description of current settings.
	"""

	network = nx.drawing.nx_agraph.to_agraph(network)
	network.draw('test.png', format='png', prog='neato')