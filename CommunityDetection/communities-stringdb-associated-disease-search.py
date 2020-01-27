import pandas as pd
import networkx as nx
import pickle
from os import path
import random
import sys

"""
This script finds all genes/proteins known to interact with a given set on input genes and the finds all HPO terms 
associated with those genes/proteins and writes it to a csv file.
"""

proG = pickle.load(open('../StringDB/protein_links_network_with_gene_names.pickle', 'rb'))
hpoG = pickle.load(open('../gene_to_hpo_with_common_names_network.pickle', 'rb'))

df = pd.DataFrame({'Community': [], 'Gene': [], 'protein_neighbor': [], 'hpo':[],
                       'hpo_name': []})
# create the blank csv with column names
df.to_csv('string-hpo-known-associations.csv')

# genes of interested
genes = ['MYLK4', 'MYH2', 'DNAH7', 'BGN', 'SAE1', 'PLK3', 'UHRF1BP1L', 'TBCB', 'UVRAG', 'ERCC3']

# communities those genes belong to based on greedy community detection
communities = [6, 6, 6, 8, 9, 9, 9, 9, 10, 10]

for i in range(len(genes)):
    gene = genes[i]
    com = communities[i]
    # given a list of genes print out all it's neighbors
    try:
        pro_neightbors = nx.neighbors(proG, gene)
    except nx.exception.NetworkXError:
        print(gene + str(' not in network'))
        continue

    out_gene = []
    out_pro_neighbor = []
    out_hpo = []
    out_name = []
    out_com = []
    for pn in pro_neightbors:

        try:
            hpo_neighbors = nx.neighbors(hpoG, pn)
        except nx.exception.NetworkXError:
            print(pn + str(' not in network'))
            continue

        for h in hpo_neighbors:
            print(i)
            out_gene.append(gene)
            out_pro_neighbor.append(pn)
            out_hpo.append(h)
            out_name.append(hpoG.nodes[h]['common_name'])
            out_com.append(com)

    # convert to a pd
    df = pd.DataFrame({'Community': out_com, 'Gene': out_gene, 'protein_neighbor': out_pro_neighbor, 'hpo': out_hpo,
                       'hpo_name': out_name})

    # append pd to csv
    df.to_csv('string-hpo-known-associations.csv', mode='a', header=False)
