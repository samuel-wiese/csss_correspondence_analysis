from networks import create_region_species_network
from plotting import draw_binary_matrix, plot_node_degrees, plot_metric_by_us_state, print_metric_by_species
from complexity import count, compute_eci, compute_fitness


geo_groups = ["United States of America"]
tax_groups = []

# Take the network
network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map, desc\
	= create_region_species_network(geo_groups, tax_groups)

# Plot obvious stuff
draw_binary_matrix(matrix_sorted, desc)
plot_node_degrees(network, desc)

# Compute metrics and draw them nicely
cou_s, cou_r = count(matrix)
eci_s, eci_r = compute_eci(matrix, regions, species)
fit_s, fit_r = compute_fitness(matrix, regions, species)
print_metric_by_species(cou_s, eci_s, fit_s, species, taxonomic_map)
plot_metric_by_us_state(cou_r, regions, desc, "Total Number of Species", "count")
plot_metric_by_us_state(eci_r, regions, desc, "ECI-Accomodatingness", "eci")
plot_metric_by_us_state(fit_r, regions, desc, "FC-Accomodatingness", "fc")
