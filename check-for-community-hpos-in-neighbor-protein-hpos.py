import pandas as pd
import networkx as nx
import pickle

"""
The purpose of this script is to check for overlap between HPO terms that appeared in community detection and HPO terms 
associated with genes/proteins that interact with the genes found in community detection. 
"""

connected_HPOs = pd.read_csv('CommunityDetection/string-hpo-known-associations.csv')
communities = pickle.load(open('CommunityDetection/greedycommunities.pickle', 'rb'))

communities_of_interest = [6, 8, 9, 10]

out_com = []
out_hpo = []
out_hpo_name = []
out_hpo_found_in = []

for i in communities_of_interest:
    com = i
    # select just the things in community i
    com_only_conencted_HPOs = connected_HPOs.loc[connected_HPOs.loc[:, 'Community'] == i, :]
    hpos = set(com_only_conencted_HPOs.loc[:, 'hpo'])
    for n in communities[i]:
        if n[0:3] == 'HP:':
            hpo = n
            # get a list of the genes/proteins that were associated with the matching HPO
            neighbors = str(
                list(com_only_conencted_HPOs.loc[com_only_conencted_HPOs.loc[:, 'hpo'] == hpo, 'protein_neighbor']))

            # get the common name of the HPO
            name = com_only_conencted_HPOs.loc[com_only_conencted_HPOs.loc[:, 'hpo'] == hpo, :]
            # check if the df is empty
            if name.shape[0] == 0:
                name = ''
            else:
                name = str(name.iloc[0, 5])

            out_com.append(com)
            out_hpo.append(hpo)
            out_hpo_name.append(name)
            out_hpo_found_in.append(neighbors)

# make a DataFrame of the accumulated lists
df = pd.DataFrame(
    {'Community': out_com, 'hpo': out_hpo, 'hpo_name': out_hpo_name, 'gene_assocaited_with_hpo': out_hpo_found_in})

# write df to a csv
df.to_csv('CommunityDetection/community-hpos-in-neighbor-protein-hpos.csv')