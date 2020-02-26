import pandas as pd
import numpy as np
from webweb import Web
import obonet
import pickle
from os import path
import networkx as nx
import copy
import seaborn as sns
import matplotlib.pyplot as plt


"""
import os
os.chdir('CommunityDetection')
"""


"""
Checks output files from the R script that ran Beckett's LPA community detection for weighted bipartite networks.
Finds the best permutation and returns communities as a 2D list.
"""


def get_communities():
    base_name = '../Becketts-raw/outfile'
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

    # fix a formatting error introduced from R
    for i in range(len(coms)):
        for j in range(len(coms[i])):
            if coms[i][j][0:3] == 'HP.':
                coms[i][j] = coms[i][j].replace('.', ':')

    return coms


"""
Returns the networkx objects of the MyGene2 data, a named and an unnamed version.
Either creates the graph from the raw tabular data or reads in a pre-made pickled version of the graph if available.
"""


def load_graphs():
    # load the MyGene2 data as a network
    g = None
    if path.exists("../../MyGene2NetworkxGraph.pickle") and path.exists("../../MyGene2NetworkxGraph-no-names.pickle"):
        print('Loading Pre-made MyGene2NetworkxGraph.pickle')
        g = pickle.load(open("../../MyGene2NetworkxGraph.pickle", "rb"))
        gn = pickle.load(open("../../MyGene2NetworkxGraph-no-names.pickle", "rb"))
    else:
        print('Creating MyGene2NetworkxGraph.pickle from stratch')
        # read in the tabular data
        data = pd.read_csv('../../Data/my_gene_2_variantes_by_family_tables.csv')

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
        pickle.dump(gn, open("../../MyGene2NetworkxGraph-no-names.pickle", "wb"))
        g = nx.relabel_nodes(g, name_mapping)
        pickle.dump(g, open("../../MyGene2NetworkxGraph.pickle", "wb"))
    return g, gn


"""
given a the MyGene2 network and the communities
plots color nodes by community and plot in WebWeb 
"""


def webweb_plot(g, coms):
    # create dictionary of node:community pairs for quick look up
    print('Plotting with WebWeb')
    coms_dict = {}
    for i in range(len(coms)):
        for c in coms[i]:
            # this is to control for string reformatting R did...
            if c[0:3] == 'HP.':
                c = c.replace('.', ':')
            coms_dict[c] = i

    # create metadata labels for plotting
    labels = {}
    for i in g.nodes:
        name = g.nodes[i]['Name']
        if name in coms_dict.keys():
            labels[i] = {'isGene': g.nodes[i]['Type'] == 'HPO', 'name': coms_dict[name], 'Community': coms_dict[name]}
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


"""
This function takes data from HPO on genes to phenotype connections, creates a networkx object of the data with the common
names of the HPO terms as metadata.
"""


def load_jenkins_gene_to_pheno():
    G = nx.DiGraph()
    if path.exists("../gene_to_hpo_with_common_names_network.pickle"):
        print('Loading Pre-made co-occurrence matrix')
        G = pickle.load(open("../gene_to_hpo_with_common_names_network.pickle", "rb"))
    else:
        # the source file was downloaded January 27th 2020 from
        # http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt\
        # Artifacts of hpo.annotations.monthly  # 159
        data = pd.read_csv('../../ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', sep='\t', skiprows=1)
        print('Generating gene_to_hpo_with_common_names_network.pickle from scratch')
        for i in range(data.shape[0]):
            # create the end
            G.add_edge(data.iloc[i, 1], data.iloc[i, 3])
            # give the hpo node a meaning full name

        nx.set_node_attributes(G, None, 'common_name')
        for i in range(data.shape[0]):
            G.nodes[data.iloc[i, 3]]['common_name'] = data.iloc[i, 2]

        pickle.dump(G, open('../gene_to_hpo_with_common_names_network.pickle', 'wb'))
    return G


"""
Take in a networkx object, return all the genes in it
"""


def get_genes(G):
    genes = []
    for n in G.nodes:
        if n[0:3] != 'HP:':
            genes.append(n)
    return genes


"""
Load the StringDB interaction network and give the node meaningful gene names, not StringDB IDs
Return the network as a networkx object
"""


def load_string_db_network():
    G = None
    if path.exists('../../StringDB/protein_interactions_named.pickle'):
        print('Loading ../../StringDB/protein_interactions_named.pickle')
        G = pickle.load(open('../../StringDB/protein_interactions_named.pickle', 'rb'))
    else:
        print('Generating ../../StringDB/protein_interactions_named.pickle\nThis may take a few minutes...')
        # create a black graph to fill
        G = nx.Graph()
        # read in the files line by line
        file = open('../../StringDB/9606.protein.actions.v11.0.txt', 'r')
        first = True
        for line in file:
            if first:
                first = False
                continue
            row = line.split('\t')
            G.add_edge(row[0], row[1])
        print('Num of nodes ' + str(len(G.nodes)))
        # load the file that contains the protein names
        pro_info = pd.read_csv('../../StringDB/protein.info.v11.0.txt', sep='\t')
        pro_info.index = pro_info.loc[:, 'protein_external_id']
        print('Number of rows: ' + str(pro_info.shape))
        print('Size of unique IDs: ' + str(len(set(pro_info.loc[:, 'protein_external_id']))))
        # create a name mapping to rename the nodes in the network
        name_mapping = {}
        for n in list(G.nodes):
            name_mapping[n] = pro_info.loc[n, 'preferred_name']
        G = nx.relabel_nodes(G, name_mapping)
        pickle.dump(G, open('../../StringDB/protein_interactions_named.pickle', 'wb'))
    return G


"""
Given the MyGene2 Gene-Phenotype graph and the communities
Return the Genes and which HPOs from their communities are not connected to them
"""


def get_genes_not_connected_to_hpos(Gj, communities):
    # this data structure will be a dictionary of dictionaries containing gene : lists of unconnected hpos
    # set up as follows:
    """
    {
        1 (a community):{
            gene: [ hpos not connected to the gene ],
            gene: [hpo1,hpo2...]
        },
        2: {
            gene: [hpo1,hpo2...],
            gene: [hpo1,hpo2...]
            ...
        }
        ...
    }
    """
    communities_genes_unconnected_hpos = {}
    com_count = -1
    not_found_count = 0
    # for each community
    for c in communities:
        com_count += 1
        # get the genes and hpos in the community
        genes = []
        hpos = []
        for n in c:
            if n[0:3] == 'HP:':
                hpos.append(n)
            else:
                genes.append(n)
        # check if a gene those hpo neighbors
        # keep track of the hpos it is not associated with
        for gene in genes:
            # get the gene's neighbors from the disease network
            try:
                neighbors = nx.neighbors(Gj, gene)
            except nx.exception.NetworkXError:
                print('Node not found: ' + str(gene))
                not_found_count += 1
                continue
            unconnected_hpos = [h for h in hpos if h not in neighbors]
            # store the unconnected HPOs only if there are some
            if len(unconnected_hpos) > 0:
                if com_count in communities_genes_unconnected_hpos.keys():
                    communities_genes_unconnected_hpos[com_count][gene] = unconnected_hpos
                else:
                    # else create a new entry in the dictionary
                    communities_genes_unconnected_hpos[com_count] = {}
                    communities_genes_unconnected_hpos[com_count][gene] = unconnected_hpos
        print('Total not found: ' + str(not_found_count))
    return communities_genes_unconnected_hpos


"""
Given the MyGene2 graph (G) and the Jenkins Gene-Phenotype graph (Gj)
Return the Genes not related to connected to HPOs in their communities by any of their neighbors in StringDB
"""


def get_genes_not_connected_to_hpos_via_neighbors(Gj, com_gene_hpo):
    print("Filtering Based on StringDB")
    # load the string DB information
    Gs = load_string_db_network()
    # for each community
    not_found_count = 0
    not_found_stringdb = 0
    found_count = 0
    # make a deep copy of the data so the keys are not being changed in the loop
    updated_com_gene_hpo = copy.deepcopy(com_gene_hpo)
    for c in com_gene_hpo.keys():
        # for each gene
        for g in com_gene_hpo[c].keys():
            # Get the gene's stringDB neighbors
            try:
                neighbors = nx.neighbors(Gs, g)
            except nx.exception.NetworkXError:
                print('Gene not found in StringDB:' + g)
                not_found_stringdb += 1
                continue
            genes_hpos = com_gene_hpo[c][g].copy()
            # for each neighbor
            for n in neighbors:
                # get the neighbor's related HPO terms
                try:
                    neighbor_hpos = nx.neighbors(Gj, n)
                except nx.exception.NetworkXError:
                    print('Gene neighbors node not found in Gj:' + n)
                    not_found_count += 1
                    continue
                else:
                    found_count += 1
                # what things are in the genes_hpos but not in neighbor_hpos? update gene_hpos
                print(list(neighbor_hpos))
                print(genes_hpos)
                print(len(genes_hpos), end=" ")

                genes_hpos = [x for x in genes_hpos if x not in list(neighbor_hpos)]
                print(len(genes_hpos))
            updated_com_gene_hpo[c][g] = genes_hpos
            # if there are not hpos left in that gene, remove it!
            if len(com_gene_hpo[c][g]) == 0:
                print('Pop goes the gene!')
                updated_com_gene_hpo[c].pop(g)
    # the number of not found may seem high, this is okay. This is the number of neighbors from StringDB not found in
    # the disease network.
    print('Total not found count: ' + str(not_found_count))
    print('Total found count: ' + str(found_count))
    print('Total not found in StringDB count: ' + str(not_found_stringdb))
    return updated_com_gene_hpo


if __name__ == "__main__":
    plt.margins(x=0.1)
    print('Running Community Detection Analysis from Beckett\'s LPAwb+ algorithm')
    # get the best best set of beckett communities
    communities = get_communities()
    # load the networkx objects a named one
    G, Gn = load_graphs()
    # plot the hair ball
    # webweb_plot(Gn, communities)
    # read in the known gene to phenotype connections
    Gj = load_jenkins_gene_to_pheno()
    # make a list of all candidate genes
    candidate_genes = get_genes(G)
    # which candidate genes are not connected to the HPOs in their community
    com_gene_hpo = get_genes_not_connected_to_hpos(Gj, communities)
    # which candidate genes are not to connected to anything their neighbors in StringDB are not connected to?
    com_gene_hpo_stringdb_filtered = get_genes_not_connected_to_hpos_via_neighbors(Gj, com_gene_hpo)

    # plot the distribution of community sizes before filtering, after step 1 of filtering and after 2 filters
    og_com_sizes = []
    filter_1_com_sizes = []
    filter_2_com_sizes = []
    # get the number of genes in each community
    for i in communities:
        len([x for x in i if x[0:3] != 'HP:'])
        og_com_sizes.append(len([x for x in i if x[0:3] != 'HP:']))
        #og_com_sizes.append(len(i))
    # get the number of genes in each community
    for i in com_gene_hpo:
        filter_1_com_sizes.append(len(com_gene_hpo[i]))

    # get the number of genes in each community
    for i in com_gene_hpo_stringdb_filtered:
        filter_2_com_sizes.append(len(com_gene_hpo_stringdb_filtered[i]))

    # Boxplot raw data
    og_name = 'Unfiltered\n N=' + str(len(og_com_sizes))
    filter_1_name = 'Direct Filtering\n N=' + str(len(filter_1_com_sizes))
    filter_2_name = 'Neighbor Filtering\n N=' + str(len(filter_2_com_sizes))
    com_lengths = og_com_sizes + filter_1_com_sizes + filter_2_com_sizes
    filter_names = [og_name] * len(og_com_sizes) + [filter_1_name] * len(filter_1_com_sizes) + [filter_2_name] * len(filter_2_com_sizes)
    data = pd.DataFrame({'Genes in each Communities': com_lengths,
                         'Filtering':filter_names })

    """
    The results seen in this figure suggest that filtering out genes based on known neighbor interactions with their 
    community's HPO terms does not reduce the number of communities. 
    """
    ax1 = sns.boxplot(x="Filtering", y="Genes in each Communities", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("filtering-communities-number-of-genes-boxplot.png")
    plt.clf()


    # plot the number of unique HPOs in each community after filtering

    og_com_sizes_g = []
    filter_1_com_sizes_g = []
    filter_2_com_sizes_g = []
    # get the number of unique HPOs in each community
    for i in communities:
        len([x for x in i if x[0:3] == 'HP:'])
        og_com_sizes_g.append(len([x for x in i if x[0:3] == 'HP:']))

    # get the number of unique HPOs in each community
    for i in com_gene_hpo:
        hpos = []
        for h in com_gene_hpo[i]:
            hpos = hpos + com_gene_hpo[i][h]
        filter_1_com_sizes_g.append(len(set(hpos)))


    # get the number of unique HPOs in each community
    for i in com_gene_hpo_stringdb_filtered:
        hpos = []
        for h in com_gene_hpo_stringdb_filtered[i]:
            hpos = hpos + com_gene_hpo_stringdb_filtered[i][h]
        filter_2_com_sizes_g.append(len(set(hpos)))


    # Boxplot of the number of unique hpos in each community
    og_name = 'Unfiltered\n N=' + str(len(og_com_sizes_g))
    filter_1_name = 'Direct Filtering\n N=' + str(len(filter_1_com_sizes_g))
    filter_2_name = 'Neighbor Filtering\n N=' + str(len(filter_2_com_sizes_g))
    com_lengths = og_com_sizes_g + filter_1_com_sizes_g + filter_2_com_sizes_g
    filter_names = [og_name] * len(og_com_sizes_g) + [filter_1_name] * len(filter_1_com_sizes_g) + [filter_2_name] * len(filter_2_com_sizes_g)
    data = pd.DataFrame({'Unique HPOs in each Communities': com_lengths,
                         'Filtering':filter_names })

    """
    The results seen in this figure suggest that filtering out genes based on known neighbor interactions with their 
    community's HPO terms does not reduce the number of communities. 
    """
    ax1 = sns.boxplot(x="Filtering", y="Unique HPOs in each Communities", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("filtering-communities-number-of-hpos-blox.png")
    plt.clf()

    ax1 = sns.violinplot(x="Filtering", y="Unique HPOs in each Communities", data=data,split=True,
                        scale="count", inner="stick")
    ax1.set_xlabel('')
    ax1.figure.savefig("filtering-communities-number-of-hpos-violin.png")
    plt.clf()

    # plot the number of unique HPOs in each GENE after filtering

    og_com_sizes_g = []
    filter_1_com_sizes_g = []
    filter_2_com_sizes_g = []
    # get the number of unique HPOs in each gene, because all genes have the same HPO at this points, I just add the legnth
    # of HPOs in the community the number of times there are genes
    for i in communities:
        for g in [x for x in i if x[0:3] != 'HP:']:
            og_com_sizes_g.append(len([x for x in i if x[0:3] == 'HP:']))

    # get the number of HPOs in each gene
    for i in com_gene_hpo:
        for h in com_gene_hpo[i]:
            hpos = com_gene_hpo[i][h]
            filter_1_com_sizes_g.append(len(set(hpos)))

    # get the number of HPOs in each gene
    for i in com_gene_hpo_stringdb_filtered:
        for h in com_gene_hpo_stringdb_filtered[i]:
            hpos = com_gene_hpo_stringdb_filtered[i][h]
            filter_2_com_sizes_g.append(len(set(hpos)))


    # Boxplot of the number of unique hpos in each community
    og_name = 'Unfiltered\n N=' + str(len(og_com_sizes_g))
    filter_1_name = 'Direct Filtering\n N=' + str(len(filter_1_com_sizes_g))
    filter_2_name = 'Neighbor Filtering\n N=' + str(len(filter_2_com_sizes_g))
    com_lengths = og_com_sizes_g + filter_1_com_sizes_g + filter_2_com_sizes_g
    filter_names = [og_name] * len(og_com_sizes_g) + [filter_1_name] * len(filter_1_com_sizes_g) + [filter_2_name] * len(filter_2_com_sizes_g)
    data = pd.DataFrame({'HPOs Associated With Each Gene': com_lengths,
                         'Filtering':filter_names })

    """
    The results seen in this figure suggest that filtering out genes based on known neighbor interactions with their 
    community's HPO terms does not reduce the number of communities. 
    """
    ax1 = sns.boxplot(x="Filtering", y="HPOs Associated With Each Gene", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("filtering-communities-hpos-per-gene-boxplot.png")
    plt.clf()

    ax1 = sns.violinplot(x="Filtering", y="HPOs Associated With Each Gene", data=data,split=True,
                        scale="count", inner="stick")
    ax1.set_xlabel('')
    ax1.figure.savefig("filtering-communities-hpos-per-gene-violin.png")
    plt.clf()
