import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from webweb import Web
import obonet
from networkx.drawing.nx_agraph import write_dot, graphviz_layout

# read in the tabular data
data = pd.read_csv('../Data/my_gene_2_variantes_by_family_tables.csv')

# get the number of genes and number of hpos
genes = list(set(data.iloc[:, 1]))
# sort them to ensure the script is deterministic as set() is not
genes.sort()

# the number of hpos and hpo names do not match. There are more names than hpo IDs
hpos = list(set(data.iloc[:, 7]))
# sort them to ensure the script is deterministic as set() is not
hpos.sort()
hpo_names = list(set(data.iloc[:, 6]))
# sort them to ensure the script is deterministic as set() is not
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
    print(i)
    G.nodes[i]['Name'] = name
    G.nodes[i]['Type'] = node_type

# find communities
c = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(G))

# create type meta data labels for webweb
labels = {}
for i in range(t):
    v = False
    n = ''
    if i < 501:
        v = True
        n = genes[i]
    else:
        n = hpos[i - 501]

    labels[i] = {'isGene': v, 'name': n}

# add the community as meta data to webweb
sub_community_count = 0
for s in c:
    for i in s:
        labels[i]['Community'] = sub_community_count
    sub_community_count += 1

w = Web(adjacency=wa.to_numpy(), display={'nodes': labels})
# set some default visualization settings in webweb
w.display.charge = 50
w.display.sizeBy = 'degree'
# genes will be colored green
w.display.colorBy = 'Community'
w.display.charge = 10
w.display.linkLength = 5
w.display.gravity = 0.5
w.show()
quit()
# webweb plot of all communities except the first way to large hair ball community
node_list = []
for i in range(1,len(c)):
    node_list = node_list + list(c[i])

H = G.subgraph(node_list)

it = 0
sublab = {}
the_map = {}
for i in list(H.nodes):
    v = False
    n = ''
    if i < 501:
        v = True
        n = genes[i]
    else:
        n = hpos[i - 501] + ' ' + data.loc[data['phenotype_hpo_id'] == hpos[i - 501]].iloc[0, 6]
    the_map[i] = it
    sublab[it] = {'isGene': v, 'name': n}
    it += 1

coms = ['zero','one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen',
        'fourteen']
# add the community as meta data to webweb
sub_community_count = 1
for s in c[1:]:
    print(s)
    print(c[sub_community_count])
    for i in s:
        print(i)
        sublab[the_map[i]]['Community'] = sub_community_count
    sub_community_count += 1
ww = Web(adjacency=nx.to_numpy_matrix(H), display={'nodes': sublab})
# set some default visualization settings in webweb
ww.display.charge = 30
ww.display.linkLength = 50
ww.display.sizeBy = 'isGene'
# genes will be colored green
ww.display.colorBy = 'Community'
ww.show()

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


shifted_G = shift_up_tree(g, G_hpo, 3, 11)
OG_shift = shifted_G.copy()
print(len(shifted_G.nodes))
c2 = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(shifted_G))
print(len(c2))
for s in list(c2):
    print(len(s))
count = 0
nodes = []
for s in list(c2):
    if count not in [0]:
        nodes = nodes + list(s)
    count+=1

# only plot nodes from the non-first communities
shifted_G = shifted_G.subgraph(nodes)

labels = {}
iter = 0
for i in shifted_G.nodes:
    v = False
    n = ''
    if shifted_G.nodes[i]['Type'] == 'HPO':
        v = False
        n = shifted_G.nodes[i]['Name']
    else:
        v = True
        n = shifted_G.nodes[i]['Name']

    com = -1
    jiter = 0
    for j in range(len(c2)):
        if i in list(c2[j]):
            com = jiter
        jiter += 1
    labels[iter] = {'isGene': v, 'name': n, 'community': com}
    iter += 1

w = Web(adjacency=nx.to_numpy_array(shifted_G), display={'nodes': labels})
# set some default visualization settings in webweb
w.display.sizeBy = 'isGene'
w.display.charge = 20
w.display.linkLength = 50
w.display.gravity = 0.05
w.display.colorBy = 'community'
w.show()


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
    return(path)


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
gg = make_network_of_tree_paths(g,G_hpo)

for n in gg.nodes:
    del gg.nodes[n]['name']

pos = graphviz_layout(gg, prog='twopi')
plt.figure(figsize=(8, 8))
nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=False)
plt.axis('equal')
plt.show()

plot = np.zeros((10000,10000))

visited = []

write_dot(gg,'test.dot')

# same layout using matplotlib with no labels
plt.title('draw_networkx')
pos =graphviz_layout(gg, prog='dot')
nx.draw(gg, pos, with_labels=False, arrows=True)
plt.show()
plt.savefig('nx_test.png')
