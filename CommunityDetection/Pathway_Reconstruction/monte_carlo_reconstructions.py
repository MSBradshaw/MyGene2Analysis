from Pathway_Reconstruction.pathway_reconstruction import count_pathways_represented_by_communities, plot_communites_and_pathways_counts
import pandas as pd
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs, \
    webweb_plot
from Expression_in_Pathways.expression_in_known_pathways import load_reactome, get_pathways
import seaborn as sns
import matplotlib.pyplot as plt
import random
random.seed(0)
# load communities
communities = get_beckett_communities()

# load reactome
reactome = load_reactome()
# get all the pathways
pathways = get_pathways(reactome, -1)
# convert pathways to dictionary of lists
pathways = {x: [y for y in pathways[x]] for x in pathways}

genes = []
for com in communities:
    genes += com
genes = [x for x in genes if x[:3] != 'HP:']

com_sizes = [len(com) for com in communities]
monte_carlo_sizes = random.choices(com_sizes,k=1000)

# make 1000 fake communities
monte_carlo_communities = []
for i in range(len(monte_carlo_sizes)):
    monte_carlo_communities.append(random.sample(genes, monte_carlo_sizes[i]))

reactome_genes = []
for p in pathways:
    reactome_genes += pathways[p]

mg2_genes_found_in_reactome = [x in reactome_genes for x in genes]
genes_names_found = [genes[i] for i in range(len(genes)) if mg2_genes_found_in_reactome[i]]

# get just the genes that are found in reactome
filtered_communities = [[gene for gene in com if gene in genes_names_found] for com in monte_carlo_communities]
# remove the empty communities
filtered_communities = [com for com in filtered_communities if len(com) > 0]

df_filtered = count_pathways_represented_by_communities(pathways, filtered_communities, False)

plot_communites_and_pathways_counts(df_filtered,'CommunityDetection/Pathway_Reconstruction/monte_carlo-filtered-communities-found-in-pathways.png')

print('Number of communities reconstructed: ' + str(len(set(df_filtered['community']))))
