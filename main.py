from networks import create_region_species_network, take_monopartite_projections
from plotting import draw_binary_matrix, plot_node_degrees, plot_metric_by_us_state, print_metric_by_species,\
	plot_clusters
from complexity import count, compute_eci, compute_fitness


# geo_groups = ["Europe"]
geo_groups = ["United-States-of-America"]
tax_groups = []

# Take the network and project it
network, matrix, matrix_sorted, regions, regions_sorted, species, species_sorted, taxonomic_map, desc\
	= create_region_species_network(geo_groups, tax_groups)
network_species, network_regions = take_monopartite_projections(network)

### Output for drawing with Gephi
nx.write_gexf(species_space, './species_space.gexf')
for cutoff in [0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 1.00]:
	MST_plus_top_links = MST_plus_links_gt_cutoff(species_space, cutoff)
	nx.write_gexf(MST_plus_top_links, f'./species_space_MST_plus_links_gt{cutoff}.gexf')

#
# # Plot obvious stuff
# draw_binary_matrix(matrix_sorted, desc)
# plot_node_degrees(network, desc)
#
# # Compute metrics and draw them nicely
# cou_s, cou_r = count(matrix)
# eci_s, eci_r = compute_eci(matrix, regions, species)
# fit_s, fit_r = compute_fitness(matrix, regions, species)
# print(regions)
# print(species)
# print(list(cou_r))
# print(list(eci_r))
# print(list(fit_r))
# print_metric_by_species(cou_s, eci_s, fit_s, species, taxonomic_map)
# plot_metric_by_us_state(cou_r, regions, desc, "Total Number of Species", "count")
# plot_metric_by_us_state(eci_r, regions, desc, "ECI-Accomodatingness", "eci")
# plot_metric_by_us_state(fit_r, regions, desc, "FC-Accomodatingness", "fc")
#
# # # Cluster species and regions
# plot_clusters(network_species, "species", desc)
# plot_clusters(network_regions, "regions", desc)

print('DONE')

# TODO: finish the similarity measures

# TODO: get a dataset on how much temperature species can take

# TODO: Analyse variance in accomodatingness and population density in mainland Europe and in China

# TODO: Build a product space, could use the backboning library
