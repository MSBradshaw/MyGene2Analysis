import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from webweb import Web
import obonet
import pickle
from os import path

"""
Checks output files from the R script that ran Beckett's LPA community detection for weighted bipartite networks.
Finds the best permutation and returns communities as a 2D list.
"""


def get_communities():
    base_name = 'CommunityDetection/Becketts/outfile'
    best_q = 0.00
    best_file = None

    # find the permutation with the best modularity score
    for i in range(1, 100):
        for j in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']:
            file = base_name + str(i) + j + '.txt'
            with open(file) as f:
                first_line = f.readline()
                if float(first_line.strip()) > best_q:
                    best_q = float(first_line.strip())
                    best_file = file

    # read in the communities associated with that best modularity score.
    # read it in as a 2D list
    first = True
    coms = []
    file = open(best_file, "r")
    for line in file:
        if first:
            first = False
        else:
            coms.append(line.strip().split(','))
    file.close()
    return coms


"""
Returns the networkx objects of the MyGene2 data, a named and an unnamed version.
Either creates the graph from the raw tabular data or reads in a pre-made pickled version of the graph if available.
"""


def load_graphs():
    # load the MyGene2 data as a network
    g = None
    if path.exists("MyGene2NetworkxGraph.pickle") and path.exists("MyGene2NetworkxGraph-no-names.pickle"):
        print('Loading Pre-made MyGene2NetworkxGraph.pickle')
        g = pickle.load(open("MyGene2NetworkxGraph.pickle", "rb"))
        gn = pickle.load(open("MyGene2NetworkxGraph-no-names.pickle", "rb"))
    else:
        print('Creating MyGene2NetworkxGraph.pickle from stratch')
        # read in the tabular data
        data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

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

        # Make networkx object
        g = nx.from_numpy_matrix(wa.to_numpy())

        # add names and types to the networkx
        nx.set_node_attributes(g, 'None', 'Name')
        nx.set_node_attributes(g, 'None', 'Type')

        name_mapping = {}
        for i in range(t):
            node_type = 'Error'
            name = 'Error'
            # upto node number 501 all nodes are genes
            if i < 501:
                node_type = 'Gene'
                name = genes[i]
            else:
                # after 501 all nodes are HPO terms
                node_type = 'HPO'
                name = hpos[i - 501]
            name_mapping[i] = name
            g.nodes[i]['Name'] = name
            g.nodes[i]['Type'] = node_type
        gn = g.copy()
        pickle.dump(gn, open("MyGene2NetworkxGraph-no-names.pickle", "wb"))
        g = nx.relabel_nodes(g, name_mapping)
        pickle.dump(g, open("MyGene2NetworkxGraph.pickle", "wb"))
    return g, gn


communities = get_communities()
G, Gn = load_graphs()


# plot them using webweb
def webweb_plot(g, coms):
    # create dictionary of node:community pairs for quick look up
    coms_dict = {}
    for i in range(len(coms)):
        for c in coms[i]:
            coms_dict[c] = i

    # create metadata labels for plotting
    labels = {}
    for i in g.nodes:
        name = g.nodes[i]['Name']
        if name in coms_dict.keys():
            labels[i] = {'isGene': g.nodes[i]['Type'] == 'HPO', 'name': name, 'Community': coms_dict[name]}
    w = Web(adjacency=nx.to_numpy_array(g), display={'nodes': labels})
    # set some default visualization settings in webweb
    w.display.charge = 50
    w.display.sizeBy = 'degree'
    # genes will be colored green
    w.display.colorBy = 'Community'
    w.display.charge = 10
    w.display.linkLength = 5
    w.display.gravity = 0.5
    w.show()


webweb_plot(Gn, communities)
# check if HPO terms in mygene2 are associated with genes or related genes
