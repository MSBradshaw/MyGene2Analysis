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

# load the pre made network of protein links
G = pickle.load( open( "protein_links_network_named.pickle", "rb" ) )

#load the protein information tsv to get meaning full names for the nodes
pro_info = pd.read_csv('StringDB/protein.info.v11.0.txt',sep='\t')

# make the row name the external id
pro_info.index = pro_info.loc[:,'protein_external_id']

g = G.copy()
# give real names to the graph
for n in G.nodes:
    g = nx.relabel_nodes(g,{n:pro_info.loc[n,'preferred_name']})
    break

pickle.dump(g,open('protein_links_networkx_gene_named.pickle','wb'))

# search for a protein/gene's associated proteins/genes
