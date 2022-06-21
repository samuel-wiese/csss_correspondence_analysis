import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap


def draw_matrix() -> None:
	"""
	Draws the matrix of species living in which regions.
	"""

	# Creates a region-species matrix
	g = nx.read_gpickle("networks/region_species.pkl")
	regions = [x for x, y in g.nodes(data=True) if y["type"] == "region"]
	species = [x for x, y in g.nodes(data=True) if y["type"] == "spicy"]
	matrix = np.zeros((len(regions), len(species)))
	for ind_reg, region in enumerate(regions):
		for ind_spe, spicy in enumerate(species):
			if [region, spicy] in g.edges:
				matrix[ind_reg, ind_spe] = 1

	# Sort it
	matrix_sorted = np.zeros((len(regions), len(species)))
	regions_ind = np.argsort(np.sum(matrix, axis=1))[::-1]
	species_ind = np.argsort(np.sum(matrix, axis=0))[::-1]
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			matrix_sorted[i, j] = matrix[regions_ind[i], species_ind[j]]

	# Draw
	plt.xlabel("Species sorted by in how many regions they are")
	plt.ylabel("Regions sorted by how\nmany species they contain")
	plt.imshow(matrix_sorted, cmap=ListedColormap(['w', 'k'], N=2), vmin=0, vmax=1, aspect=1)
	plt.tight_layout()
	plt.savefig("img/matrix.pdf", dpi=300, bbox_inches="tight")
	plt.clf()


draw_matrix()
