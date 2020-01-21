import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from webweb import Web
import obonet
import pickle
from os import path
from math import sqrt
import random

random.seed(0)

def calc_phi(tp, fp, tn, fn):
    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)


# read in the tabular data
data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

# load the pre-made and labeled weighted bipartite of the MyGene2 data
G = pickle.load(open('MyGene2NetworkxGraph.pickle', "rb"))

# do greedy community detection
communities = nx.algorithms.community.modularity_max.greedy_modularity_communities(G)

# g will be the network with combine community nodes
g = G.copy()
# for each of the interesting communities
for i in [6, 8, 9, 10]:
    # get all the nodes in the community that are HPOs
    hpos = []
    for n in communities[i]:
        if G.nodes[n]['Type'] == 'HPO':
            hpos.append(n)
    # merge all HPOs to be one node
    g = nx.contracted_nodes(g, hpos[0], hpos[1])
    g.nodes[hpos[0]]['Name'] = 'Community-' + str(i) + g.nodes[hpos[0]]['Name'] + ',' + hpos[1]
    for j in range(2, len(hpos)):
        g = nx.contracted_nodes(g, hpos[0], hpos[j])
        g.nodes[hpos[0]]['Name'] = g.nodes[hpos[0]]['Name'] + ',' + hpos[j]
    print(g.nodes[hpos[0]]['Name'])

A = nx.to_numpy_array(g)

for e in g.edges:
    if g.get_edge_data(e[0], e[1])['weight'] > 1:
        print(e)
        print(g.get_edge_data(e[0], e[1])['weight'])
        print('greater')
        break

# count the number of genes in the network
genes = []
hpos = []
for n in g.nodes:
    if g.nodes[n]['Type'] == 'Gene':
        genes.append(n)
    else:
        hpos.append(n)

# create a co occurance matrix from the adjacency matrix, rows are genes, cols and HPOs
co_matrix = A[0:len(genes), len(genes):A.shape[1]]

co = pd.DataFrame(co_matrix)

co.columns = hpos
co.index = genes
com_phis = []
for i in [6, 8, 9, 10]:
    print('Community: ' + str(i))
    # get the HPOs and the first gene
    com_hpo = None
    com_gene = []
    for n in list(communities[i]):
        if G.nodes[n]['Type'] != 'HPO':
            com_gene.append(n)
        elif com_hpo is None:
            # if it is the first occurrence of an HPO save it
            com_hpo = n
    # calc phi for the community and each gene
    for g in com_gene:
        true_pos = co.loc[g, com_hpo]
        false_neg = sum(co.loc[g, :]) - true_pos
        false_pos = sum(co.loc[:, com_hpo]) - true_pos
        true_neg = sum(sum(co_matrix)) - sum(co.loc[g, :]) - sum(co.loc[:, com_hpo]) + co.loc[g, com_hpo]
        p = calc_phi(true_pos, false_pos, true_neg, false_neg)
        com_phis.append(p)
        print(com_hpo + ' ' + g + ' : ' + str(p))


def combine_nodes_phi(graph, gene, combine_hpos):
    # merge all HPOs to be one node
    graph = nx.contracted_nodes(graph, combine_hpos[0], combine_hpos[1])
    graph.nodes[combine_hpos[0]]['Name'] = 'Community-' + str(i) + graph.nodes[combine_hpos[0]]['Name'] + ',' + \
                                           combine_hpos[1]
    for j in range(2, len(combine_hpos)):
        graph = nx.contracted_nodes(graph, combine_hpos[0], combine_hpos[j])
        graph.nodes[combine_hpos[0]]['Name'] = graph.nodes[combine_hpos[0]]['Name'] + ',' + combine_hpos[j]

    # create adjacency matrix
    Adj = nx.to_numpy_array(graph)

    # count the number of genes in the network and store them for use as col and row names
    genes = []
    hpos = []
    for n in graph.nodes:
        if graph.nodes[n]['Type'] == 'Gene':
            genes.append(n)
        else:
            hpos.append(n)

    # create a co occurrence matrix from the adjacency matrix, rows are genes, cols and HPOs
    co_matrix = Adj[0:len(genes), len(genes):Adj.shape[1]]
    co = pd.DataFrame(co_matrix)
    co.columns = hpos
    co.index = genes

    # calc phi
    com_hpo = combine_hpos[0]
    true_pos = co.loc[gene, com_hpo]
    false_neg = sum(co.loc[gene, :]) - true_pos
    false_pos = sum(co.loc[:, com_hpo]) - true_pos
    true_neg = sum(sum(co_matrix)) - sum(co.loc[gene, :]) - sum(co.loc[:, com_hpo]) + co.loc[gene, com_hpo]
    p = calc_phi(true_pos, false_pos, true_neg, false_neg)
    return p


# now boot strap (if it is fair to call it that...)
# select a random gene with multiple neighbors
# combine up to 6 (max number of HPOs in interesting communities) of it's neighbors into a super node
# calc phi for this random gene and it's super neighbor
# see if community 8's phi falls above the 95% threshold

phis = []
for i in range(10):
    gene = genes[random.randint(0, len(genes))]
    # get neighbors of gene, because it is bipartite
    neighbors = nx.neighbors(G, gene)
    neighbors = list(neighbors)
    # randomly choose how many neighbors will be combine (3-6)
    number_to_choose = random.randint(3, 7)
    if number_to_choose > len(neighbors):
        # there are not enough neighbors
        print('skipping')
        continue
    # choose the indexes of neighbors to be combine
    neigh_indexes = random.sample(range(0, len(list(neighbors))), number_to_choose)

    # extact the names of those neighbors
    chosen_neighbors = []
    for j in neigh_indexes:
        chosen_neighbors.append(neighbors[j])

    # calc phi
    phis.append(combine_nodes_phi(G.copy(), gene, chosen_neighbors))
print(phis)

# TODO 1000+ boot straps and plot the data. Box and Manhattan Plots
