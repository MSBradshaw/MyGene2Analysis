from os import path
import pickle
import networkx as nx
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import random
import copy
import time
import multiprocessing as mp
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_communities, load_graphs, \
    load_string_db_network, load_jenkins_gene_to_pheno
from sklearn.decomposition import PCA

# just temporarily
GTEX_GLOBAL = None

"""
import os
os.chdir('CommunityDetection/Expression-in-Pathways')
"""

"""
Load the Reactome Data with gene common names mapping to a reactome pathway
"""


def load_reactome():
    reactome_path = "../../Reactome/reactome.pickle"
    reactome = None
    if path.exists(reactome_path):
        print('Loading reactome')
        reactome = pickle.load(open(reactome_path, 'rb'))
    else:
        print('Generating Reactome from scratch')
        # make a dictionary of uniprot ids and gene names
        uniprot_dict = {}
        with open("../../Uniprot/HUMAN_9606_id_gene_name_mapping.txt") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                row = line.split(sep='\t')
                uniprot_dict[row[0]] = row[2].strip()
        reactome = nx.Graph()
        key_erros_count = 0
        with open("../../Reactome/UniProt2Reactome_All_Levels_Human_Only.txt") as fp:
            while True:
                line = fp.readline()
                if not line:
                    break
                row = line.strip().split(sep='\t')
                # add an edge from the uniprot id (mapped to the common gene name) and the pathway
                try:
                    unip_id = row[0]
                    if '-' in unip_id:
                        unip_id = re.sub('-\d+', '', unip_id)
                    reactome.add_edge(uniprot_dict[unip_id], row[1])
                    reactome.nodes[uniprot_dict[unip_id]]['Type'] = 'Gene'
                    reactome.nodes[uniprot_dict[unip_id]]['Info'] = 'Gene'
                    reactome.nodes[row[1]]['Type'] = 'Pathway'
                    reactome.nodes[row[1]]['Info'] = row[3]
                except KeyError:
                    print('Key Error: ' + row[0])
                    key_erros_count += 1
        print(reactome.nodes['R-HSA-109582'])
        print('Key Error Count: ' + str(key_erros_count))

    return reactome


"""
Returns a pandas DataFrame of the GTEx gene expression data
"""


def load_gtex():
    gtex_pickle_path = '../../GTEx/GTEx.pickle'
    gtex = None
    if path.exists(gtex_pickle_path):
        print('loading GTEx')
        gtex = pickle.load(open(gtex_pickle_path, 'rb'))
    else:
        print('Reading GTEx: This could take a while...')
        gtex = pd.read_csv('../../GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=2)
        # protocol 4 is to allow for large file handling (over 4gb)
        pickle.dump(gtex, open(gtex_pickle_path, 'wb'), protocol=4)
    return gtex


"""
Return a 2D list of pathways from reactome
"""


def get_pathways(reactome, number_of_pathways, seed=0):
    # get all pathways
    pathways = {}
    for node in reactome.nodes:
        if reactome.nodes[node]['Type'] == 'Pathway':
            genes = nx.neighbors(reactome, node)
            pathways[node] = genes
    if number_of_pathways == -1:
        print('returning all pathways')
        return pathways
    # choose 20 pathways at "random"
    random.seed(seed)
    # select random key indexes
    randos = [random.randint(0, len(pathways)) for p in range(0, number_of_pathways)]
    path_keys = list(pathways.keys())
    # get the keys at those indexes
    selected_keys = [path_keys[r] for r in randos]
    # get the genes associated with each key (pathway)
    selected_paths = {}
    for key in selected_keys:
        selected_paths[key] = list(pathways[key])
    return selected_paths


"""
Given the reactome networkx object
Create and save 2 box plots of the sizes (number of genes) of the pathways
"""


def plot_pathway_sizes(reactome):
    sizes = []
    for node in reactome.nodes:
        if reactome.nodes[node]['Type'] == 'Pathway':
            # get the number of neighbors
            neighbors = list(nx.neighbors(reactome, node))
            num_of_neighs_that_are_pathways = len([n for n in neighbors if reactome.nodes[n]['Type'] == 'Pathway'])
            if num_of_neighs_that_are_pathways > 0:
                # this will print out the number of pathways associated with each pathway,
                # it never prints so we know the network is a true bipartite
                print(num_of_neighs_that_are_pathways)
            sizes.append(len(neighbors))
    print(sizes)
    data = pd.DataFrame({'Number of Genes': sizes, 'All Pathways': ['All Pathways'] * len(sizes)})

    # box plot
    ax1 = sns.boxplot(x="All Pathways", y="Number of Genes", data=data)
    ax1.set_xlabel('')
    ax1.figure.savefig("size-of-reactome-pathways.png")
    plt.clf()

    # box plot log scale
    ax1 = sns.boxplot(x="All Pathways", y="Number of Genes", data=data)
    ax1.set_xlabel('')
    ax1.set_yscale("log")
    ax1.figure.savefig("size-of-reactome-pathways-log.png")
    plt.clf()
    print('Mean:' + str(np.mean(np.array(sizes))))
    print(len(sizes))


"""
Given a pandas DataFrame
Return the same DataFrame but without the rows labeled delete
"""


def remove_filler_rows(df):
    return df[df['pathway'] != 'delete']


"""
Given a dictionary of pathways and the GTEx data 
Return a Pandas DataFrame formated like so:
Gene,   Pathway,    Tissue1,    Tissue2, ...
SEMA4A, R-HSA-109582,   867,    5309, ...
UBA6,   R-HSA-109582,   442,    0446, ...
"""


def get_tabular_expression_data(pathways):
    file = ''
    if isinstance(pathways, list):
        file = pathways[1]
        pathways = pathways[0]
        if path.exists(file):
            print('File Exists: ' + file)
            return
    else:
        print('not list' + str(type(pathways)))

    gtex = GTEX_GLOBAL
    gtex.index = gtex['Description']
    # get make of copy of gtex to copy all of the columns but get 2 rows (and remove those rows later)
    exp_data = copy.deepcopy(gtex.iloc[0:2, :])
    exp_data.insert(1, 'pathway', ['delete', 'delete'])
    not_in_gtex_count = 0
    for p in pathways.keys():
        print('P: ' + str(p))
        for gene in pathways[p]:
            # get each of its expression values
            try:
                try:
                    gene_exp = gtex.loc[gene]
                    gene_exp.loc['pathway'] = p
                    exp_data.loc[exp_data.shape[0]] = gene_exp
                except KeyError:
                    print(gene + " not found in GTEx")
                    not_in_gtex_count += 1
            except ValueError:
                print('Row Error: ' + str(p))
                return None

    print('Not found in GTEx Count: ' + str(not_in_gtex_count))
    print("pickling " + file)
    pickle.dump(remove_filler_rows(exp_data), open(file, 'wb'), protocol=4)
    return None


"""
Given the reactome networkx object and the gtex Pandas DataFrame
Return a Pandas DataFrame of the Gene Expression data of all Genes, this is a large file, over 4gb
"""


def get_expression_for_all_genes(reactome):
    # -1 will return all pathways
    paths = get_pathways(reactome, -1)
    expression_pickle_path = 'all-genes-expression-pd.pickle'
    if path.exists(expression_pickle_path):
        print('Loading Gene Expression From Pickle')
        return pickle.load(open(expression_pickle_path, 'rb'))
    else:
        print('Generating Gene Expression From Scratch this may take a while...')
        start = time.time()
        # do this in parallel, with all but 2 threads available (so I can still do other things while it runs)
        pool = mp.Pool(processes=40)
        # create pairs of parameters to be used for each parallel task,
        # each one will get 1 pathway as a dictionary and True
        # (which will tell the function to use the global variable EXPRESSION_GLOBAL instead of returning)
        param_lists = [[{key: paths[key]}, key + ".pickle"] for key in paths]
        # run jobs in parallel
        pool.map(get_tabular_expression_data, param_lists)
        files = [key + ".pickle" for key in paths]
        exp = None
        for f in files:
            if exp is None:
                exp = pickle.load(open(f, 'rb'))
            else:
                try:
                    exp = pd.concat([exp, pickle.load(open(f, 'rb'))], ignore_index=True, sort=False)
                except FileNotFoundError:
                    print('File Does Not Exist: ' + f)
        # remove the dumpy filler data
        exp = exp.loc[exp['pathway'] != 'delete']
        print('done')
        # global EXPRESSION_GLOBAL
        print(exp.shape)
        print('Generating took ' + str(time.time() - start))
        # protocol 4 is so pickle can file objects over the size of 4gb
        exp = remove_filler_rows(exp)
        pickle.dump(exp, open(expression_pickle_path, 'wb'), protocol=4)
        return exp


"""
Given the reactome
Return the median gene expressions per tissue table from GTEx but of only the genes found in MyGene2
Deletes from memory GTEX_GLOBAL
"""


def get_reduced_gtex(reactome):
    filepath = 'gtex-pandas-mygene2-genes-only.pickle'
    global GTEX_GLOBAL
    if path.exists(filepath) and False:
        print('Loading Reduced GTEX')
        GTEX_GLOBAL = None
        return pickle.load(open(filepath, 'rb'))
    else:
        print('Generating Reduced GTEX')
        exp = get_expression_for_all_genes(reactome)
        # load mygene2
        G, Gn = load_graphs()
        genes = [x for x in G.nodes if x[0:3] != 'HP:']
        exp.index = exp['Description']
        GTEX_GLOBAL = None
        exp = exp.loc[genes]
        exp = remove_filler_rows(exp)
        pickle.dump(exp, open(filepath, 'wb'), protocol=4)
        return exp


"""
Given a pandas DataFrame with samples as rows and variables as columns
Return a pandas DataFrame of the first 5 principle components
"""


def get_principle_components(df):
    pca = PCA(n_components=5)
    df = df.fillna(0)
    data = df.iloc[:, 3:].T
    pca.fit(data)
    pc_df = pd.DataFrame(pca.components_.T)
    pc_df.columns = ['PC1 ' + str(round(pca.explained_variance_ratio_[0] * 100, 2)) + '%',
                     'PC2 ' + str(round(pca.explained_variance_ratio_[1] * 100, 2)) + '%',
                     'PC3 ' + str(round(pca.explained_variance_ratio_[2] * 100, 2)) + '%',
                     'PC4 ' + str(round(pca.explained_variance_ratio_[3] * 100, 2)) + '%',
                     'PC5 ' + str(round(pca.explained_variance_ratio_[4] * 100, 2)) + '%']
    pc_df.insert(0, 'Name', df.index)
    pc_df.insert(1, 'pathway', df['pathway'].to_list())
    pc_df.insert(2, 'Description', df['Description'].to_list())
    return pc_df


"""
Given the reactome 
Returns a pandas DataFrame of the first 5 principle components from the reduced gtex data
"""


def get_reduced_gtex_pca(reactome):
    filepath = 'gtex-reduced-pandas-pca.pickle'
    if path.exists(filepath) and False:
        print('Loading Reduced PCA')
        return pickle.load(open(filepath, 'rb'))
    else:
        print('Creating Reduced PCA')
        exp = get_reduced_gtex(reactome)
        e2 = get_principle_components(exp)
        pickle.dump(e2, open(filepath, 'wb'), protocol=4)
        return e2


"""
Given the reactome
Returns a pandas DataFrame of the first 5 principle components from the gtex data
"""


def get_gtex_pca(reactome):
    filepath = 'gtex-pandas-pca.pickle'
    if path.exists(filepath) and False:
        print('Loading PCA')
        return pickle.load(open(filepath, 'rb'))
    else:
        print('Creating PCA')
        exp = get_expression_for_all_genes(reactome)
        pca = PCA(n_components=5)
        exp = exp.fillna(0)
        data = exp.iloc[:, 3:].T
        pca.fit(data)
        e2 = pd.DataFrame(pca.components_.T)
        e2.columns = ['PC1 ' + str(round(pca.explained_variance_ratio_[0] * 100, 2)) + '%',
                      'PC2 ' + str(round(pca.explained_variance_ratio_[1] * 100, 2)) + '%',
                      'PC3 ' + str(round(pca.explained_variance_ratio_[2] * 100, 2)) + '%',
                      'PC4 ' + str(round(pca.explained_variance_ratio_[3] * 100, 2)) + '%',
                      'PC5 ' + str(round(pca.explained_variance_ratio_[4] * 100, 2)) + '%']
        e2.insert(0, 'Name', exp.index)
        e2.insert(1, 'pathway', exp['pathway'].to_list())
        e2.insert(2, 'Description', exp['Description'].to_list())
        pickle.dump(e2, open(filepath, 'wb'), protocol=4)
        return e2


"""
Given a set of pathways and the reactome networkx object
Make a PCA plot of the given expression data
"""


def pca_plot(pathways, reactome):
    print('Making PCA Plot')
    # load each individual community in
    raw_gtex = None
    for p in pathways:
        if raw_gtex is None:
            raw_gtex = pickle.load(open('Pathway-Pickles/' + p + '.pickle', 'rb'))
        else:
            raw_gtex = raw_gtex.append(pickle.load(open('Pathway-Pickles/' + p + '.pickle', 'rb')))
    raw_gtex = remove_filler_rows(raw_gtex)
    # get the PC's of the raw_gtex data
    pca = get_principle_components(raw_gtex)
    pca.insert(3, 'pathway-name', 'No Info')
    # get the reactome info for each path way so they have meaningful names
    # create a mapping of the names
    pathway_names = {}
    for p in set(pca['pathway']):
        pathway_names[p] = reactome.nodes[p]['Info']
    # assign the meaning full names
    print(pca.columns)
    for i in range(pca.shape[0]):
        print(pca.iloc[i, 1])
        pca.iloc[i, 3] = pathway_names[pca.iloc[i, 1]]
    # the data from pca that is path of the pathways of interest
    data = pca[pca['pathway'] == pathways[0]]
    for i in range(1, len(pathways)):
        data = data.append(pca[pca['pathway'] == pathways[i]])
    # plot it!
    ax = sns.scatterplot(x=data.columns[4], y=data.columns[5], data=data, hue="pathway-name")
    plt.show()
    plt.clf()
    ax = sns.scatterplot(x=data.columns[4], y=data.columns[6], data=data, hue="pathway-name")
    plt.show()


if __name__ == "__main__":
    reactome = load_reactome()
    GTEX_GLOBAL = load_gtex()
    # plot_pathway_sizes(reactome)
    paths = get_pathways(reactome, 5, 2)
    path_names = list(paths.keys())
    # e = get_reduced_gtex(reactome)
    # pc = get_reduced_gtex_pca(reactome)
    pca_plot(path_names, reactome)
