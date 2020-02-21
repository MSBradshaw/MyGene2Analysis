from becketts_community_analysis import get_communities, load_graphs
from webweb import Web
import networkx as nx
import numpy as np

"""
import os
os.chdir('CommunityDetection')
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
    w = Web(adjacency=nx.to_numpy_array(subgraph),display={'nodes': labels})
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
Given the named networkx object and the community of interest (list)
Returns the "strong" one-mode projection of the network
When I say strong I mean that all genes in the network as neighbors of each other in a one-mode projection
"""


def strong_one_mode_projection(G,community):
    # make a dict of each gene and its neighbors
    gene_neighbors = {}



if __name__ == "__main__":
    print('Properties of Communities')
communities = get_communities()
G, Gn = load_graphs()

plot_communities_one_by_one(Gn, communities)

plot_single_community_webweb(Gn,communities,10)

get_gene_neighbors_of_gene(G, 'SATB2', communities[10])
