import pandas as pd
import networkx as nx
import pickle
from os import path
import random
import sys
# load the pre made network of protein links
G = pickle.load( open( "protein_links_network_named.pickle", "rb" ) )

#load the protein information tsv to get meaning full names for the nodes
pro_info = pd.read_csv('StringDB/protein.info.v11.0.txt',sep='\t')

# make the row name the external id
pro_info.index = pro_info.loc[:,'protein_external_id']

g = G.copy()

name_mapping = {}
print(sys.argv)
start = int(sys.argv[1]) * 1000
end = start + 1000
# give real names to the graph
for n in list(G.nodes)[start:end]:
    print(n)
    name_mapping[n] = pro_info.loc[n,'preferred_name']

# g = nx.relabel_nodes(g,{n:pro_info.loc[n,'preferred_name']})
pickle.dump(name_mapping,open('name_mapping_'+str(start/1000)+'.pickle','wb'))
# pickle.dump(g,open('protein_links_networkx_gene_named.pickle','wb'))

# search for a protein/gene's associated proteins/genes
