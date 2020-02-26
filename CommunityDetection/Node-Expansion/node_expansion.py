from Becketts_Community_Analysis_Results.becketts_community_analysis import get_communities, load_graphs, \
    load_string_db_network, load_jenkins_gene_to_pheno
import networkx as nx
import numpy as np
import copy
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

"""
import os
os.chdir('CommunityDetection/Node-Expansion')
"""

"""
Given a network of genes and hpo terms, return a list of the number of HPOs associated with each gene
"""


def get_hpos_per_gene(G):
    hpos_per_gene = []
    for n in G.nodes:
        if n[0:3] != 'HP:':
            hpos_per_gene.append(len(list(nx.neighbors(G, n))))
    return hpos_per_gene


"""
Given a network of genes and hpo terms
Returns a list of the number of HPOs associated with each gene and it's neighbors in StringDB and 
a list containing Number of Genes from given network found in StringDB, # not found in string DB and 
number of gene neighbors from StringDB found in the given network and number not found
"""


def get_hpos_per_expanded_gene(G, name, verbose=True):
    Gs = load_string_db_network()
    hpos_per_gene = []
    not_found_neighbor_count = 0
    found_neighbor_count = 0
    given_graph_not_found_count = 0
    given_graph_found_count = 0
    for n in G.nodes:
        if n[0:3] != 'HP:':
            # get it's neighbors in StringDB
            try:
                neighbors = nx.neighbors(Gs, n)
            except nx.exception.NetworkXError:
                if verbose:
                    print(name + ' Gene not found in StringDB: ' + n)
                given_graph_not_found_count += 1
            else:
                given_graph_found_count += 1
            # get the HPOs associated with each neighbor and the gene itself
            hpos = list(nx.neighbors(G, n))
            for neighbor in neighbors:
                try:
                    hpos = hpos + list(nx.neighbors(G, neighbor))
                except nx.exception.NetworkXError:
                    if verbose:
                        print('Neighbor not in ' + name + ': ' + '|' + neighbor + '|')
                    not_found_neighbor_count += 1
                else:
                    found_neighbor_count += 1

            # convert the list to a set to get only unique hpos, add number of unique hpos to list
            hpos_per_gene.append(len(set(hpos)))
    print('\n')
    print(name + ' Found in StringDB: ' + str(given_graph_found_count))
    print(name + ' Not Found in StringDB: ' + str(given_graph_not_found_count))
    # the vast majority of StringDB neighbors are not found in MyGene2
    print('Neighbors Found: ' + str(found_neighbor_count))
    print('Neighbors Not Found: ' + str(not_found_neighbor_count))
    return hpos_per_gene, [given_graph_found_count,
                           given_graph_not_found_count,
                           found_neighbor_count,
                           not_found_neighbor_count]


if __name__ == "__main__":
    print('Node Expansion Analysis')
    G, Gn = load_graphs()
    Gj = load_jenkins_gene_to_pheno()
    # get the number of HPOs associated with each MyGene2 gene
    mg2_hpos_counts = get_hpos_per_gene(G)
    print('MyGene2 average number of HPOs per gene: ' + str(np.mean(np.array(mg2_hpos_counts))))
    # get the number of HPOs associated with each Jenkins gene
    jenkins_hpos_counts = get_hpos_per_gene(Gj)
    print('Jenkins average number of HPOs per gene: ' + str(np.mean(np.array(jenkins_hpos_counts))))
    # get the number of HPOs associated with each MyGene2 gene after expansion with StringDB data
    mg2_expanded_hpo_counts, mg2_stats = get_hpos_per_expanded_gene(G, 'MyGene2', False)
    print('MyGene2 Expanded average number of HPOs per gene: ' + str(np.mean(np.array(mg2_expanded_hpo_counts))))
    # get the number of HPOs associated with each expanded Jenkins gene
    jenkins_expanded_hpos_counts, jenkins_stats = get_hpos_per_expanded_gene(Gj, 'Jenkins', False)
    print('Jenkins Expanded average number of HPOs per gene: ' + str(np.mean(np.array(jenkins_expanded_hpos_counts))))

    # produce data
    counts = mg2_hpos_counts + jenkins_hpos_counts + mg2_expanded_hpo_counts + jenkins_expanded_hpos_counts
    groups = ['MyGene2'] * len(mg2_hpos_counts) + \
             ['Jenkins'] * len(jenkins_hpos_counts) + \
             ['MyGene2\nExpanded'] * len(mg2_expanded_hpo_counts) + \
             ['Jenkins\nExpanded'] * len(jenkins_expanded_hpos_counts)

    data = pd.DataFrame({'HPOs per Gene': counts, 'Group': groups})

    # box plot
    ax1 = sns.boxplot(x="Group", y="HPOs per Gene", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("hpos-per-gene-boxplot.png")
    plt.clf()

    # box plot with log scale
    ax1 = sns.boxplot(x="Group", y="HPOs per Gene", data=data)
    ax1.set_xlabel('')
    ax1.set_yscale("log")
    ax1.figure.savefig("log-hpos-per-gene-boxplot.png")
    plt.clf()

    # histogram
    ax1 = sns.boxplot(x="Group", y="HPOs per Gene", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("hpos-per-gene-boxplot.png")
    plt.clf()

    # histogram log scale
    ax1 = sns.barplot(x="Group", y="HPOs per Gene", data=data)
    ax1.set_xlabel('')
    ax1.set_yscale("log")
    ax1.figure.savefig("log-hpos-per-gene-histogram.png")
    plt.clf()

    # plot the stats about the expansion
    counts = mg2_stats[0:2] + jenkins_stats[0:2]
    groups = ['MyGene2', 'MyGene2', 'Jenkins', 'Jenkins']
    types = ['Found', 'Not Found', 'Found', 'Not Found']
    data = pd.DataFrame({'Group': groups, 'Count': counts, 'Type': types})

    ax1 = sns.barplot(x="Group", y="Count", hue='Type', data=data)
    ax1.set_xlabel('')
    ax1.set_ylabel('Genes Not Found In StringDB')
    # ax1.set_yscale("log")
    ax1.figure.savefig("genes-not-found-in-StringDB-histogram.png")
    plt.clf()

    # plot the stats about the expansion and number of neighbors not found in MyGene2 and Jenkins
    counts = mg2_stats[2:] + jenkins_stats[2:]
    groups = ['MyGene2', 'MyGene2', 'Jenkins', 'Jenkins']
    types = ['Found', 'Not Found', 'Found', 'Not Found']
    data = pd.DataFrame({'Group': groups, 'Count': counts, 'Type': types})


    ax1 = sns.barplot(x="Group", y="Count", hue='Type', data=data)
    ax1.set_xlabel('')
    ax1.set_ylabel('StringDB Neighbors Not Found')
    ax1.figure.subplots_adjust(left=0.2)
    ax1.figure.savefig("neighbors-not-found-histogram.png")
    plt.clf()
