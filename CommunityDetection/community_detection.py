import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from webweb import Web
import obonet
import pickle
from os import path

# read in the tabular data
data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

# get the number of genes and number of hpos and hpo name (the later 2 do not match, FYI)
genes = list(set(data.iloc[:, 1]))
hpos = list(set(data.iloc[:, 7]))
hpo_names = list(set(data.iloc[:, 6]))

# sort them to ensure the script is deterministic as set() is not
genes.sort()
hpos.sort()
hpo_names.sort()

# make a weighed adjacency matrix starting from an all zeros pandas
t = len(genes) + len(hpos)
wa = pd.DataFrame(np.zeros((t, t)), columns=genes + hpos, index=genes + hpos)

# add weight to all connected edges
for i in range(data.shape[0]):
    row = data.iloc[i, :]
    wa.loc[row[1], row[7]] += 1

# MODULARITY ANALYSIS
# Make networkx object
G = nx.from_numpy_matrix(wa.to_numpy())

# add names and types to the networkx
nx.set_node_attributes(G, 'None', 'Name')
nx.set_node_attributes(G, 'None', 'Type')

for i in range(t):
    node_type = 'Error'
    name = 'Error'
    if i < 501:
        node_type = 'Gene'
        name = genes[i]
    else:
        node_type = 'HPO'
        # if it is an HPO term then set the name equal to the term id and semantic name
        name = hpos[i - 501]  # + ' ' + data.loc[data['phenotype_hpo_id'] == hpos[i - 501]].iloc[0, 6]
    G.nodes[i]['Name'] = name
    G.nodes[i]['Type'] = node_type

# read in the hpo obo file
url = '/Users/michael/PycharmProjects/MyGene2Analysis/hp.obo'
g = obonet.read_obo(url)

dists = []
for h in hpos:
    current = h
    dist = 0
    try:
        while g.nodes[current]['name'] != 'All':
            current = g.nodes[current]['is_a'][0]
            dist += 1
        dists.append(dist)
    except KeyError:
        print('Key Error' + current + ' dist: ' + str(dist))

# make a box plot of the distances
plt.figure()
plt.boxplot(dists, 1)
plt.savefig('Figures/HPO-distance-box-plot.png')

# the attribute is_a move up the tree toward more general terms
def move_up_tree(graph, hpo, move_up, min_dist):
    c = hpo
    d = 0
    path = [hpo]
    try:
        while graph.nodes[c]['name'] != 'All':
            c = graph.nodes[c]['is_a'][0]
            d += 1
            path.append(c)
    except KeyError:
        print('Key Error ' + current + ' dist: ' + str(dist))
        return None
    if min_dist == None or d > min_dist:
        if len(path) - move_up > min_dist:
            return path[move_up]
        else:
            # return the item at min_dist
            print('Cannot move ' + str(move_up) + ' spaces up the tree. Returning node at minimum distance')
            return path[-(min_dist + 1)]
    else:
        print('Path shorter than min distance. Returning starting point')
        return hpo


def shift_up_tree(hpo_graph, mg2_graph, move_up, min_dist):
    # Generate a dictionary of node that are suppose to be combine
    new_nodes_dict = {}
    for h in list(mg2_graph.nodes):
        if mg2_graph.nodes[h]['Type'] != 'HPO':
            continue
        else:
            hpo = mg2_graph.nodes[h]['Name']
        new_node = move_up_tree(hpo_graph, hpo, move_up, min_dist)
        if new_node in new_nodes_dict.keys():
            new_nodes_dict[new_node].append(hpo)
        else:
            new_nodes_dict[new_node] = [hpo]
    new_graph = mg2_graph.copy()
    # loop over dict and combind all nodes
    for k in new_nodes_dict.keys():
        if k is None:
            print('k is None')
            continue
        if len(new_nodes_dict[k]) == 1:
            # just rename the one node
            new_graph = nx.relabel_nodes(new_graph, {new_nodes_dict[k][0]: k + '-new'})
        else:
            # if we get in here we assume that the number of nodes to be renamed is >= 2
            # if rename the first node to be the new name
            combine_name = new_graph.nodes[new_nodes_dict[k][0]]['Name'] + ',' + new_graph.nodes[new_nodes_dict[k][1]][
                'Name']
            new_graph = nx.relabel_nodes(new_graph, {new_nodes_dict[k][0]: k + '-new'})
            # combind the first two nodes
            new_graph = nx.contracted_nodes(new_graph, k + '-new', new_nodes_dict[k][1])
            new_graph.nodes[k + '-new']['Name'] = combine_name
            # iteratively combind all other nodes in the list
            for i in range(2, len(new_nodes_dict[k])):
                combine_name = new_graph.nodes[k + '-new']['Name'] + ',' + new_graph.nodes[new_nodes_dict[k][i]]['Name']
                new_graph = nx.contracted_nodes(new_graph, k + '-new', new_nodes_dict[k][i])
                new_graph.nodes[k + '-new']['Name'] = combine_name
    # now change all of the names to not include -new
    for n in new_graph.nodes:
        if new_graph.nodes[n]['Type'] == 'HPO':
            new_graph = nx.relabel_nodes(new_graph, {n: n.replace('-new', '')})
    return new_graph


# build a version of G where the node is actually labels with the HBO term
G_hpo = G.copy()
for n in G_hpo.nodes:
    G_hpo = nx.relabel_nodes(G_hpo, {n: G_hpo.nodes[n]['Name']})

# save the weighted and named bipartite for later use in other scripts
pickle.dump(G_hpo, open('MyGene2NetworkxGraph.pickle','wb'))

communities = nx.algorithms.community.modularity_max.greedy_modularity_communities(G_hpo)
nodes =  ['CNOT1', 'AGT', 'HP:0000707', 'HP:0001085', 'ANK2', 'ATXN1', 'DBH', 'HP:0001658', 'GALT', 'SIRPG', 'MCF2L2', 'HP:0002046', 'HP:0002027', 'HP:0100633', 'HP:0001880', 'MEGF9', 'HP:0003002', 'HP:0000975', 'HP:0000131', 'PAPOLG', 'HP:0002664', 'TSFM', 'NCOA1', 'HP:0002094', 'HP:0010536', 'HP:0012432', 'HP:0002315', 'HP:0001336', 'SPR', 'HP:0002758', 'HP:0000006', 'HP:0011134', 'HP:0100013', 'DALRD3', 'CDK5RAP3', 'HP:0001382', 'ATRAID', 'MICAL2', 'AGRN', 'NAGA', 'ATXN', 'SHKBP1', 'FOXE1', 'HP:0002893', 'HP:0100702', 'ITPKC', 'FAM186A', 'HP:0001909', 'HP:0001029', 'HP:0012532', 'HP:0002270', 'HP:0001271', 'ANP32B', 'F10', 'DACH1', 'HP:0000718', 'PIK3C2A', 'HP:0001369', 'HP:0010783', 'HP:0000007', 'HP:0001955', 'HP:0000618', 'MUC16', 'HP:0007021', 'MSRA', 'SLC12A6', 'HP:0003737', 'HP:0012514', 'PPP1R26', 'DMXL1', 'HP:0100790', 'SEC14L5', 'AACS', 'HCLS1', 'KIF1B', 'ZNRF4', 'KSR2', 'GABRG2', 'HP:0001539', 'C5ORF30', 'MORN3', 'HP:0000668', 'HP:0000766', 'HP:0000023', 'HP:0003270', 'HP:0002558', 'ZIC1', 'TCN2', 'DNAJC13', 'SNAPC1', 'HP:0100272', 'UBA2', 'KLF6', 'HP:0000411', 'FREM2', 'SPTB', 'ZNF227', 'HP:0000776', 'BCAM', 'ANKRD46', 'TMEM108', 'SCARA3', 'PKD1', 'KCNH2', 'HP:0001697', 'ACTG2', 'ATAD3B', 'ZNF534', 'HP:0000089', 'HP:0100771', 'FBXO11', 'SCN1A', 'ARMC3', 'HP:0025004', 'HP:0002511', 'HP:0001863', 'HP:0002656', 'HP:0001025', 'HP:0025230', 'HP:0012317', 'HP:0002654', 'HP:0000722', 'HP:0009926', 'HP:0030883', 'HP:0010886', 'HP:0030881', 'ACAN', 'HP:0001370', 'HP:0002652', 'HP:0001864', 'HP:0011971', 'HP:0001944', 'HP:0040292', 'HP:0007565', 'DNER', 'HP:0001560', 'HP:0001787', 'HP:0012194', 'SHANK2', 'HP:0010780', 'HP:0100033', 'HP:0011443', 'HP:0002031', 'HP:0000736', 'HP:0002037', 'HP:0006530', 'AGR2', 'HP:0000121', 'HP:0006685', 'LRP1', 'HP:0012387', 'HP:0004469', 'HP:0005567', 'HP:0001644', 'HP:0001989', 'HP:0000282', 'MYLK4', 'TMEM39B', 'MYH2', 'DNAH7', 'HP:0002803', 'HP:0002202', 'HP:0000546', 'HP:0007675', 'HP:0000608', 'HP:0000512', 'CEP250', 'HP:0001133', 'HP:0007924', 'HP:0000456', 'HP:0001829', 'HP:0002984', 'BGN', 'HP:0010503', 'HP:0010442', 'HP:0009556', 'SAE1', 'HP:0002665', 'HP:0002671', 'PLK3', 'UHRF1BP1L', 'TBCB', 'HP:0030430', 'HP:0000613', 'HP:0000563', 'HP:0000992', 'UVRAG', 'HP:0030127', 'ERCC3', 'ZMYM2', 'HP:0030736', 'HP:0009792', 'HP:0011787', 'HP:0000826', 'HP:0001973', 'HP:0000872', 'BCL6B', 'HP:0001596', 'HP:0012377', 'HP:0020046', 'MAP1S', 'HP:0030516']

degs = []
for n in nodes:
    degs.append(G_hpo.degree[n])

counts = []
g_list = list(data.iloc[:, 1])
h_list = list(data.iloc[:, 7])
for n in nodes:
    if n[0:3] == 'HP:':
        counts.append(h_list.count(n))
    else:
        counts.append(g_list.count(n))

num_coms_shift = []
if path.exists("num_coms_shift.pickle1"):
    print('Loading Pre-made num_coms_shift')
    num_coms_shift = pickle.load(open("num_coms_shift.pickle", "rb"))
else:
    # how does changing the shift effect the number of communities? Hold min distance constant at 6 (the average distance)
    for i in range(20,21):
        print('Shift: ' + str(i))
        shifted_G = shift_up_tree(g, G_hpo, i, 6)
        c = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(shifted_G))
        num_coms_shift.append(len(c))

    pickle.dump(num_coms_shift,open('num_coms_shift.pickle','wb'))

num_coms_min_dist = []
if path.exists("num_coms_min_dist.pickle1"):
    print('Loading Pre-made num_coms_min_dist')
    num_coms_min_dist = pickle.load(open("num_coms_min_dist.pickle", "rb"))
else:
# how does changing the minimum distance effect the number of communities? Hold shift constant at 3 (arbitrary)
    for i in range(20,21):
        print('Min Dist: ' + str(i))
        shifted_G = shift_up_tree(g, G_hpo, 3, i)
        c = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(shifted_G))
        num_coms_min_dist.append(len(c))

    pickle.dump(num_coms_min_dist,open('num_coms_min_dist.pickle','wb'))


x = range(1,21)
plt.plot(x, num_coms_min_dist, label="Minimum Distance")
plt.plot(x, num_coms_shift, label="Shift Distance")
plt.legend()

fig = plt.figure()
ax = fig.gca()
ax.set_xticks(np.arange(1, 21, 1))
ax.set_yticks(np.arange(1, 21, 1))
ax.set_xlabel('Shift Distance / Minimum Distance')
ax.set_ylabel('Number of Communities')
plt.plot(x, num_coms_min_dist, label="Minimum Distance")
plt.plot(x, num_coms_shift, label="Shift Distance")
plt.legend()
plt.grid()
plt.savefig('Figures/community-detection-shifts-min-distance-effect.png')


