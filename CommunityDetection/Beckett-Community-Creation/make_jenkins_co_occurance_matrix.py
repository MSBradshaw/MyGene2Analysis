from Becketts_Community_Analysis_Results.becketts_community_analysis import load_jenkins_gene_to_pheno
import networkx as nx
import pandas as pd

# load Jenkins
jenkins = load_jenkins_gene_to_pheno()

# create a co-occurrence matrix with HPO columns and gene rows

# get all genes and all hpos from Jenkins
genes = []
hpos = []
for n in list(jenkins.nodes):
    if n[:3] == 'HP:':
        hpos.append(n)
    else:
        genes.append(n)

co_dict = {h:[] for h in hpos}

# build the matrix row by row
for g in genes:
    neighbors = list(nx.neighbors(jenkins,g))
    for h in hpos:
        # if the genes if neighbors with the hpo term, append 1 else 0
        if h in neighbors:
            co_dict[h].append(1)
        else:
            co_dict[h].append(0)

# conver to a pandas DataFrame
df = pd.DataFrame(co_dict)
# create row indexes
df.index = genes
# save to csv for use in R or else where
df.to_csv('CommunityDetection/Beckett-Community-Creation/jenkins-co_occurrence_matrix.csv')