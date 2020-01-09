import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from webweb import Web
import obonet
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import pygraphviz

# read in the tabular data
data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')
# get the number of genes and number of hpos
genes = list(set(data.iloc[:, 1]))
# the number of hpos and hpo names do not match. There are more names than hpo IDs
hpos = list(set(data.iloc[:, 7]))
hpo_names = list(set(data.iloc[:, 6]))
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
            print(current)
            dist += 1
        dists.append(dist)
        print(dist)
    except KeyError:
        print('Key Error' + current + ' dist: ' + str(dist))


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
    print(path)
    if min_dist == None or d > min_dist:
        print(len(path))
        if len(path) - move_up > min_dist:
            return path[move_up]
        else:
            # return the item at min_dist
            print('Cannot move ' + str(move_up) + ' spaces up the tree. Returning node at minimum distance')
            return path[-(min_dist + 1)]
    else:
        print('Path shorter than min distance. Returning starting point')
        return hpo


# build a version of G where the node is actually labels with the HBO term
G_hpo = G.copy()
for n in G_hpo.nodes:
    G_hpo = nx.relabel_nodes(G_hpo, {n: G_hpo.nodes[n]['Name']})


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


def tree_path(graph, hpo):
    print(hpo)
    c = hpo
    d = 0
    path = [hpo]
    try:
        while graph.nodes[c]['name'] != 'All':
            c = graph.nodes[c]['is_a'][0]
            d += 1
            path.append(c)
            print(c)
    except KeyError:
        print('Key Error ' + current + ' dist: ' + str(dist))
        return None
    return (path)


def make_network_of_tree_paths(hpo_graph, mg2_graph):
    paths = []
    for node in mg2_graph.nodes:
        if mg2_graph.nodes[node]['Type'] == 'HPO':
            new_paths = tree_path(hpo_graph, node)
            if new_paths is not None:
                paths = new_paths + paths
    list(set(paths))
    return hpo_graph.subgraph(list(set(paths)))


# Generate just the portion of the HPO tree that we are using
gg = make_network_of_tree_paths(g, G_hpo)


def traverse(graph, dictionary, current, depth):
    paths = graph.in_edges(current)
    accessible = len(paths)
    for p in paths:
        dictionary, a = traverse(graph, dictionary, p[0], depth + 1)
        accessible += a
    dictionary[depth].append(accessible)
    return dictionary, accessible


dicts = {}

for i in range(1000):
    dicts[i] = []

dicts, a = traverse(gg, dicts, 'HP:0000001', 1)

rows = []
for k in dicts.keys():
    rows.append()
    if len(dicts[k]) > 0:
        print([i for i in dicts[k] if i != 0 and i != 1])