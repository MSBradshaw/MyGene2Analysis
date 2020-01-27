import networkx as nx
import pickle

G = pickle.load( open( "protein_links_network_named.pickle", "rb" ) )

name_mapping = {}

# load name mappings broken up into multiple files
# these files were generated by XXXXX.py
for i in range(0,21):
    name_mapping.update(pickle.load(open('name_mapping_' + str(i) + '.0.pickle','rb')))

# rename the nodes
G = nx.relabel_nodes(G,name_mapping)

# save the named stringDB network
pickle.dump(G,open('StringDB/protein_links_network_with_gene_names.pickle','wb'))