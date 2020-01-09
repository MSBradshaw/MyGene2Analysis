import networkx as nx
import pandas as pd
import pickle

data = pd.read_csv('StringDB/9606.protein.links.v11.0.txt', sep=' ')
name_data = pd.read_csv('StringDB/protein.info.v11.0.txt', sep='\t')

G = nx.Graph()

for i in range(data.shape[0]):
    row = data.iloc[i, :]
    # the weight is the combined_score from the file, not sure what is means...
    G.add_edge(row[0], row[1], weight=int(row[2]))

# read in the normal edge list

pickle.dump(G, open("protein_links_network.pickle", "wb"))

G = pickle.load( open( "protein_links_network.pickle", "rb" ) )

nx.set_node_attributes(G, '', 'Name')
success = 0
fail = 0
for i in range(name_data.shape[0]):
    row = name_data.iloc[i, :]
    identifier = row[0]
    gene_name = str(row[1]).upper()
    try:
        try:
            G.nodes[identifier]['Name'] = gene_name
            success += 1
        except KeyError:
            fail += 1
    except AttributeError:
        fail += 1

pickle.dump(G, open("protein_links_network_named.pickle", "wb"))
