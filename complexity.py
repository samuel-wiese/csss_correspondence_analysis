import numpy as np

from typing import Tuple, List, Optional


def count(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
	"""
	Simply counts the number of regions each species is not in (pickyness), and the number of species each region
	contains.

	Parameters
	----------
	matrix : np.ndarray
		Matrix of species and regions.

	Returns
	-------
	cou_s, cou_r : Tuple[np.ndarray, np.ndarray]
		Simply count.
	"""

	return matrix.shape[0] - matrix.sum(axis=0), matrix.sum(axis=1)


def compute_eci(matrix: np.ndarray, regions: List[str], species: List[str]) -> Tuple[np.ndarray, np.ndarray]:
	"""
	Computes the HH ECI.

	Parameters
	----------
	matrix : np.ndarray
		Matrix of species and regions.
	regions : List[str]
		The list of regions.
	species : List[str]
		The list of species.

	Returns
	-------
	eci_s, eci_r : Tuple[np.ndarray, np.ndarray]
		The ECI for species and for regions.
	"""

	# Compute the ECI for species and for regions
	pickyness = matrix.shape[0] - matrix.sum(axis=0)
	pickyness[np.where(pickyness == 0)] = 1e-12
	accomodatingness = matrix.sum(axis=1)

	# Calculate ECI and PCI eigenvectors
	mcp1 = matrix / accomodatingness[:, np.newaxis]
	mcp2 = matrix / pickyness[np.newaxis, :]
	mcp2_t = mcp2.T.copy()
	matrix_r = mcp1 @ mcp2_t
	matrix_s = mcp2_t @ mcp1
	r_w, r_v = np.linalg.eig(matrix_r)
	s_w, s_v = np.linalg.eig(matrix_s)
	s_v = np.real(s_v[:, s_w.argsort()[::-1]])[:, 1]
	r_v = np.real(r_v[:, r_w.argsort()[::-1]])[:, 1]
	s_v = (s_v - np.mean(s_v)) / np.std(s_v)
	r_v = (r_v - np.mean(r_v)) / np.std(r_v)

	# Fix signs
	if "Florida" in regions:
		if r_v[regions.index("Florida")] < 0:
			r_v = -r_v
	if "Egretta" in species:
		if s_v[species.index("Egretta")] > 0:
			s_v = -s_v
	if "Alcedo" in species:
		if s_v[species.index("Alcedo")] > 0:
			s_v = -s_v
	if "Anas" in species:
		if s_v[species.index("Anas")] > 0:
			s_v = -s_v

	return s_v, r_v


def compute_fitness(matrix: np.ndarray, regions: List[str], species: List[str], extremality: Optional[float] = 1.0)\
		-> Tuple[np.ndarray, np.ndarray]:
	"""
	Computes the Pietronero fitness.

	Parameters
	----------
	matrix : np.ndarray
		Matrix of species and regions.
	regions : List[str]
		The list of regions.
	species : List[str]
		The list of species.
	extremality : Optional[float] = 1.0
		Extremality parameters governing the influence of the least accommodating country.

	Returns
	-------
	fit_s, fit_r : Tuple[np.ndarray, np.ndarray]
		The FC for species and for regions.
	"""

	# Compute the fitness of regions (accomodatingness) and the complexity (pickyness) of species
	pickyness = np.ones(matrix.shape[1])
	accomodatingness = np.ones(matrix.shape[0])
	tol = 1
	while tol > 0.0001:
		new_acc = matrix.dot(pickyness)
		new_acc /= new_acc.mean()
		new_pickyness = (1 / matrix.T.dot(1 / accomodatingness ** extremality)) ** (- 1 / extremality)
		new_pickyness /= new_pickyness.mean()
		tol = (np.abs(pickyness - new_pickyness).mean() + np.abs(accomodatingness - new_acc).mean()) / 2
		pickyness = new_pickyness.copy()
		accomodatingness = new_acc.copy()

	# Normalize them
	accomodatingness = (accomodatingness - np.mean(accomodatingness)) / np.std(accomodatingness)
	pickyness = (pickyness - np.mean(pickyness)) / np.std(pickyness)

	# Fix signs
	if "Florida" in regions:
		if accomodatingness[regions.index("Florida")] < 0:
			accomodatingness = -accomodatingness
	if "Egretta" in species:
		if pickyness[species.index("Egretta")] > 0:
			pickyness = -pickyness
	if "Alcedo" in species:
		if pickyness[species.index("Alcedo")] > 0:
			pickyness = -pickyness
	if "Anas" in species:
		if pickyness[species.index("Anas")] > 0:
			pickyness = -pickyness

	return pickyness, accomodatingness
