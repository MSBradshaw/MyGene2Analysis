from Becketts_Community_Analysis_Results.becketts_community_analysis import get_communities, load_graphs
from webweb import Web
import networkx as nx
import numpy as np
import copy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

"""
import os
os.chdir('CommunityDetection/Strong-One-Mode-Filtering-Results')
"""

"""
Given the networkx object (with nodes numbered, not named), the community 2D list and the index of the community of 
interest
"""


def plot_single_community_webweb(g, coms, com):
    print('Plotting with WebWeb')

    # create dictionary of node:community pairs for quick look up
    coms_dict = {}
    for i in range(len(coms)):
        for c in coms[i]:
            coms_dict[c] = i

    # get the indexes of the nodes in the community of interest
    com_indexes = []
    for i in g.nodes:
        name = g.nodes[i]['Name']
        if name in coms_dict.keys():
            if coms_dict[name] == com:
                com_indexes.append(i)

    # get a subset of the grapha that is just the community of interest
    subgraph = g.subgraph(com_indexes)

    # create a name mapping to relabel subgraph nodes starting at one
    name_mapping = {}
    node_count = 0
    for i in subgraph.nodes:
        print(i)
        name_mapping[i] = node_count
        node_count += 1

    # relabel subgraph nodes starting at one
    subgraph = nx.relabel_nodes(subgraph, name_mapping)

    # create labels for the nodes for use in webweb
    labels = {}
    for i in subgraph.nodes:
        name = subgraph.nodes[i]['Name']
        if name in coms_dict.keys():
            if coms_dict[name] == com:
                labels[i] = {'isGene': subgraph.nodes[i]['Type'] == 'HPO', 'name': name,
                             'Community': coms_dict[name]}

    # plot it using webweb
    w = Web(adjacency=nx.to_numpy_array(subgraph), display={'nodes': labels})
    # set some default visualization settings in webweb
    w.display.charge = 50
    w.display.sizeBy = 'degree'
    w.display.colorBy = 'isGene'
    w.display.charge = 20
    w.display.linkLength = 50
    w.display.gravity = 0.1
    w.show()


"""
Given the networkx object (with nodes numbered, not named) and the community 2D list
Will go through each community plotting it in webweb but wait for user input (even just an enter key) before plotting 
the next one
"""


def plot_communities_one_by_one(Gn, coms):
    for i in range(len(communities)):
        plot_single_community_webweb(Gn, communities, i)
        g = input("Press Anything to continue: ")


"""
Given the named networkx object, the gene of interest and the gene's community
returns a list of all genes in the community that are 1 hop neighbor sof the gene of interest
"""


def get_gene_neighbors_of_gene(G, gene, community):
    # get all the HPO terms associated with the gene (because it is a bipartite)
    hpos = nx.neighbors(G, gene)
    hpos = list(hpos)
    # get all the genes connected to each of it's HPOs
    genes = []
    for h in hpos:
        neighbors = list(nx.neighbors(G, h))
        genes = genes + neighbors
    # get just the unique elements from the genes list
    genes = np.unique(np.array(genes)).tolist()
    # return just the 1 hop neighbors that are in the community
    return [x for x in genes if x in community]


"""
Given a dictionary, sorting it by the length of teh values
Taken from: https://stackoverflow.com/questions/16868457/python-sorting-dictionary-by-length-of-values
But it has been modified...
"""


def sort_by_values_len(dict):
    dict_len = {key: len(value) for key, value in dict.items()}
    import operator
    sorted_key_list = sorted(dict_len.items(), key=operator.itemgetter(1), reverse=True)
    sorted_dict = {item[0]: dict[item[0]] for item in sorted_key_list}
    return sorted_dict


"""
Given the named networkx object and the community of interest (list)
Returns the "strong" one-mode projection of the network
When I say strong I mean that all genes in the network are neighbors of each other in a one-mode projection
"""


def strong_one_mode_projection(G, community):
    # make a dict of each gene and its neighbors
    gene_neighbors = {}
    for node in community:
        if node[0:3] != 'HP:':
            # it is a gene
            gene_neighbors[node] = get_gene_neighbors_of_gene(G, node, community)
    # sort the dictionary
    gene_neighbors = sort_by_values_len(gene_neighbors)
    changed = True
    com_genes = [x for x in community if x[0:3] != 'HP:']
    # while there are still nodes being remove from the matrix
    while changed:
        changed = False
        new_gene_neighbors = copy.deepcopy(gene_neighbors)
        # for each node in the dictionary
        for node in gene_neighbors:
            # get the nodes that are in the community but not in the gene's neighborhood
            not_in_neighborhood = [x for x in com_genes if x not in gene_neighbors[node]]
            # destroy nodes not in the neighborhood
            for gene in not_in_neighborhood:
                if gene in new_gene_neighbors.keys():
                    new_gene_neighbors.pop(gene)
                    changed = True
        gene_neighbors = copy.deepcopy(new_gene_neighbors)

    # re-put together the community but with on the real strong genes
    hpos = [x for x in community if x[0:3] == 'HP:']
    new_com = hpos + list(gene_neighbors.keys())

    # make a subgraph of the community
    sub = nx.subgraph(G,new_com)
    # filter out HPOs that do not connect genes
    hpos = [x for x in hpos if sub.degree(x) > 1]

    # make final community
    new_com = hpos + list(gene_neighbors.keys())

    return new_com


"""
Given the named networkx object and 2D community list
Return a filtered version of communities where only genes that are all neighbors or each other remain
"""


def strong_one_mode_filter(G, communities):
    new_communities = copy.deepcopy(communities)
    for i in range(len(communities)):
        com = communities[i]
        print("Com: " + str(i))
        print('Pre: ' + str(len([x for x in com if x[0:3] != 'HP:'])))
        new_communities[i] = strong_one_mode_projection(G, com)
        print('Post: ' + str(len([x for x in new_communities[i] if x[0:3] != 'HP:'])))
    return new_communities


"""
Given a 2D list of the communities pre and post filtering to a strong 1-mode projection
Makes 2 violin plots of the information, one of genes and one of HPOs
"""


def plot_communities_sizes_pre_and_post_filtering(communities_pre, communities_post):
    print('plotting community sizes')
    # get the number of genes in each community but filter out communites with less than 2 genes in them
    pre_sizes = [len([node for node in c if node[0:3] != 'HP:']) for c in communities_pre
                 if len([node for node in c if node[0:3] != 'HP:']) >= 2]
    post_sizes = [len([node for node in c if node[0:3] != 'HP:']) for c in communities_post
                  if len([node for node in c if node[0:3] != 'HP:']) >= 2]

    # get the number of HPOs in each community
    pre_sizes_hpo = [len([node for node in c if node[0:3] == 'HP:']) for c in communities_pre
                 if len([node for node in c if node[0:3] == 'HP:']) >= 2]
    post_sizes_hpo = [len([node for node in c if node[0:3] == 'HP:']) for c in communities_post
                  if len([node for node in c if node[0:3] == 'HP:']) >= 2]

    print(pre_sizes)
    print(post_sizes)

    pre_name = 'Pre Filtering\n N=' + str(len(pre_sizes))
    post_name = 'Post Filtering\n N=' + str(len(post_sizes))

    com_lengths = pre_sizes + post_sizes
    filter_names = [pre_name] * len(pre_sizes) + [post_name] * len(post_sizes)
    data = pd.DataFrame({'Genes in each Communities': com_lengths,
                         'Filtering': filter_names})

    ax1 = sns.violinplot(x="Filtering", y="Genes in each Communities", data=data, split=True,
                         scale="count", inner="stick")
    ax1.set_xlabel('')
    ax1.figure.savefig("strong-one-mode-filter-Gene-violin.png")
    plt.clf()


    pre_name = 'Pre Filtering\n N=' + str(len(pre_sizes_hpo))
    post_name = 'Post Filtering\n N=' + str(len(post_sizes_hpo))

    com_lengths = pre_sizes_hpo + post_sizes_hpo
    filter_names = [pre_name] * len(pre_sizes_hpo) + [post_name] * len(post_sizes_hpo)
    data = pd.DataFrame({'HPO in each Communities': com_lengths,
                         'Filtering': filter_names})

    ax1 = sns.violinplot(x="Filtering", y="HPO in each Communities", data=data, split=True,
                         scale="count", inner="stick")
    ax1.set_xlabel('')
    ax1.figure.savefig("strong-one-mode-filter-HPO-violin.png")
    plt.clf()


if __name__ == "__main__":
    print('Properties of Communities')
    communities = get_communities()
    G, Gn = load_graphs()

    # plot_communities_one_by_one(Gn, communities)
    c2 = strong_one_mode_filter(G, communities)
    plot_single_community_webweb(G,communities,10)
    plot_single_community_webweb(G, c2, 10)
    plot_communities_sizes_pre_and_post_filtering(communities, c2)
