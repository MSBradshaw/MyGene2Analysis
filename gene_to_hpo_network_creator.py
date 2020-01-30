import pandas as pd
import networkx as nx
import pickle
from os import path
import random
import sys

"""
This file takes data from HPO on genes to phenotype connections, creates a networkx object of the data with the common
names of the HPO terms as metadata.
"""

# the source file was downloaded January 27th 2020 from
# http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt

data = pd.read_csv('ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', sep='\t', skiprows=1)

G = nx.DiGraph()

for i in range(data.shape[0]):
    # create the end
    G.add_edge(data.iloc[i, 1], data.iloc[i, 3])
    # give the hpo node a meaning full name

nx.set_node_attributes(G, None, 'common_name')
for i in range(data.shape[0]):
    G.nodes[data.iloc[i, 3]]['common_name'] = data.iloc[i, 2]

pickle.dump(G, open('gene_to_hpo_with_common_names_network.pickle', 'wb'))
