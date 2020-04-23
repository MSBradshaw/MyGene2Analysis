import pandas as pd
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs, \
    webweb_plot
from Expression_in_Pathways.expression_in_known_pathways import load_reactome, get_pathways
import seaborn as sns
import matplotlib.pyplot as plt

# load communities
communities = get_beckett_communities()

# load reactome
reactome = load_reactome()
# get all the pathways
pathways = get_pathways(reactome, -1)
# convert pathways to dictionary of lists
pathways = {x: [y for y in pathways[x]] for x in pathways}

"""
Given:
pathways - a dictionary of list representing a pathway and the genes in it
community - list of genes and hpo terms (hpo terms are not used, but are controlled for)
strict_threshold - if True reconstructions only count if all genes are present in the pathway
strict_threshold - if False reconstructions count if half or more genes are present in a pathway
Return a list of pathways that the community is found in (found meaning less than half the genes are missing)
"""


def find_pathways_matching_communities(pathways, community, strict_threshold=False):
    found_pathways = []
    com = [x for x in community if x[:3] != 'HP:']
    threshold = 1
    if not strict_threshold:
        threshold = int(len(community) / 2)
    for p in pathways:
        missing_count = 0
        for gene in com:
            if gene not in pathways[p]:
                missing_count += 1
        if missing_count < threshold:
            found_pathways.append(p)
            # print(pathways[p])
            # print(missing_count)
    return found_pathways


"""
Given:
pathways - a dictionary of list representing a pathway and the genes in it
coms - list lists of genes and hpo terms (hpo terms are not used, but are controlled for)
strict_threshold - if True reconstructions only count if all genes are present in the pathway
strict_threshold - if False reconstructions count if half or more genes are present in a pathway
Return pandas DataFrame with the following columns ['community_size', 'pathway_size', 'community', 'pathway']
"""


def count_pathways_represented_by_communities(pathways, coms, strict_threshold):
    size_of_com = []
    size_of_path = []
    com_names = []
    path_names = []
    communities_reconstructed_count = 0
    for i in range(len(coms)):
        paths = find_pathways_matching_communities(pathways, coms[i], strict_threshold)
        if len(paths) > 0:
            for p in paths:
                size_of_com.append(len(filtered_communities[i]))
                com_names.append(str(i) + ' ')
                size_of_path.append(len(pathways[p]))
                path_names.append(p)
            communities_reconstructed_count += len(paths)
            print(str(i) + ' community size: ' + str(len(filtered_communities[i])) + ' ' + str(paths))
    print('Sum of pathways communities are present in: ' + str(communities_reconstructed_count))
    return pd.DataFrame({'community_size': size_of_com,
                         'pathway_size': size_of_path,
                         'community': com_names,
                         'pathway': path_names})


def plot_communites_and_pathways_counts(df,name):
    count = 0
    for c in list(set(df['community'])):
        colors = ['r', 'g', 'b', 'orange', 'purple', 'black', 'yellow']
        sub = df[df['community'] == c]
        plt.scatter(sub['community_size'], sub['pathway_size'], color=colors[count])
        count += 1
    plt.xlabel('Community Size')
    plt.ylabel('Pathways Size')
    plt.savefig(name)


com_df = count_pathways_represented_by_communities(pathways, communities, False)
plot_communites_and_pathways_counts(com_df,'CommunityDetection/Pathway-Reconstruction/communities-found-in-pathways.png')

# how many of the MyGene2 Genes are in reactome
genes = []
for com in communities:
    genes += com
genes = [x for x in genes if x[:3] != 'HP:']

reactome_genes = []
for p in pathways:
    reactome_genes += pathways[p]

mg2_genes_found_in_reactome = [x in reactome_genes for x in genes]
print('MyGene2 genes found in reactome: ' + str(sum(mg2_genes_found_in_reactome)))
mg2_genes_not_found_in_reactome = [x not in reactome_genes for x in genes]
print('MyGene2 genes Not found in reactome: ' + str(sum(mg2_genes_not_found_in_reactome)))

genes_names_not_found = [genes[i] for i in range(len(genes)) if not mg2_genes_found_in_reactome[i]]
genes_names_found = [genes[i] for i in range(len(genes)) if mg2_genes_found_in_reactome[i]]

# get just the genes that are found in reactome
filtered_communities = [[gene for gene in com if gene in genes_names_found] for com in communities]
# remove the now empty communites
filtered_communities = [com for com in filtered_communities if len(com) > 0]

df_filtered = count_pathways_represented_by_communities(pathways, filtered_communities, False)

plot_communites_and_pathways_counts(df_filtered,'CommunityDetection/Pathway-Reconstruction/filtered-communities-found-in-pathways.png')

df_filtered_strict = count_pathways_represented_by_communities(pathways, filtered_communities, True)
plot_communites_and_pathways_counts(df_filtered_strict,'CommunityDetection/Pathway-Reconstruction/strict_filtered-communities-found-in-pathways.png')