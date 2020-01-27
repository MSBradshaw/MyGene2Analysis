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


name_mapping = {}

# the first argument of the script should be which node to start on
# creating the mapping is broken up into chunks of 1000 so it can be run in parallel
# this script will need to be run for inputs [0-20] to generate all mappings
start = int(sys.argv[1]) * 1000
# the script should create mappings for 1000 nodes at a time
end = start + 1000

# give real names to the graph
for n in list(G.nodes)[start:end]:
    print(n)
    name_mapping[n] = pro_info.loc[n,'preferred_name']

#save this portion of the remapping names
pickle.dump(name_mapping,open('name_mapping_'+str(start/1000)+'.pickle','wb'))
