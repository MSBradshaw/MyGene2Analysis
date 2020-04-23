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
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs
import seaborn as sns

random.seed(0)

def calc_phi(tp, fp, tn, fn):
    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)


# read in the tabular data
data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

# load the weighted bipartite of the MyGene2 data
G, Gn = load_graphs()

# do greedy community detection
communities = get_beckett_communities()

# g will be the network with combine community nodes
g = G.copy()
# for each of the interesting communities
for i in range(len(communities)):
    # get all the nodes in the community that are HPOs
    hpos = []
    for n in communities[i]:
        if G.nodes[n]['Type'] == 'HPO':
            hpos.append(n)
    # check there are more than 1 HPO, if not there is no need to merge nodes
    if len(hpos) == 1:
        continue
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
for i in range(len(communities)):
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
gene_phis = {}
while len(phis) < 10000:
    if len(phis) % 100 == 0:
        print(len(phis))
    try:
        index = random.randint(0, len(genes)-1)
        gene = genes[index]
    except IndexError:
        print(index)
        print(len(genes))
        quit()
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
    p = combine_nodes_phi(G.copy(), gene, chosen_neighbors)
    phis.append(p)
    if gene in gene_phis.keys():
        gene_phis[gene].append(p)
    else:
        gene_phis[gene] = [p]

print(len(phis))
pickle.dump(phis, open('CommunityDetection/10000-bootstraps-phi-becketts.pickle','wb'))
pickle.dump(gene_phis, open('CommunityDetection/10000-bootstraps-phi-dictionary-becketts.pickle','wb'))
pickle.dump(com_phis,open('CommunityDetection/beckett-community-phis.pickle','wb'))

phis = pickle.load(open('CommunityDetection/10000-bootstraps-phi-becketts.pickle', 'rb'))
gene_phis = pickle.load(open('10000-bootstraps-phi-dictionary-becketts.pickle', 'rb'))
com_phis = pickle.load(open('CommunityDetection/beckett-community-phis.pickle', 'rb'))

df = pd.DataFrame({'phi':phis + com_phis,'type':(['Boot Strap']*len(phis)) + (['Communities'] * len(com_phis)) })

sns.boxplot(y="phi",x='type',data=df)
plt.savefig('CommunityDetection/phi-becketts-bootstrap-10000-boxplots.png')

#
# fig2, ax2 = plt.subplots()
# # plot the only community with an intersting phu (com 8)
# # plt.axhline(y=0.31613179478426595,color='blue',linestyle=':')
# # plt.text(1, 0.31613179478426595, '     Community 8', fontsize=8, va='bottom', ha='left', backgroundcolor=(0.0, 0.0, 0.0, 0.0))
# # plot the threshold for significance using a Bon Ferroni correction (we did 11 tests of communities)
# # thresh = np.percentile(np.array(phis),(100 - (5 / 11)))
#
# # plt.axhline(y=thresh,color='red',linestyle='--')
# # plt.text(1,thresh, '     99.54 Threshold', fontsize=8, va='bottom', ha='left', backgroundcolor=(0.0, 0.0, 0.0, 0.0))
# ax2.set_title('Phi\'s Bootstraping, n=10,000')
# ax2.boxplot(phis, notch=True)
# fig2.savefig('CommunityDetection/phi-10000-boxplot-becketts.png', dpi=fig2.dpi)
