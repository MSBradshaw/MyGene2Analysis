from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs, \
    load_string_db_network, load_jenkins_gene_to_pheno
import pandas as pd

beckett_communities = get_beckett_communities()

tsg = pd.read_csv('TissueSpecificGenes/tissue_specific_genes.csv')


# how many tissues are there that are specific to every gene in the community
def get_tissues_specific_to_all_genes(tissue_df,com):
    # get just the genes in the community
    genes = [x for x in com if x[:3] != 'HP:']
    if len(genes) <= 1:
        # print('Warning, only one gene in the community')
        return []
    tissues = []
    # find the tissues that gene is specific for
    tissue_specific_genes_count = 0
    for g in genes:
        temp_tissues = list(tissue_df[tissue_df['Description'] == g]['tissue'])
        if len(temp_tissues) > 0:
            tissue_specific_genes_count += 1
        tissues += temp_tissues

    if tissue_specific_genes_count == 1:
        # print('Only one tissues specific gene')
        return []

    common_tissues = []
    for t in list(set(tissues)):
        if tissues.count(t) >= tissue_specific_genes_count:
            common_tissues.append(t)
    return common_tissues


# how many tissues are there that are specific to every gene in the community
def get_tissues_specific_to_half_genes(tissue_df,com):
    # get just the genes in the community
    genes = [x for x in com if x[:3] != 'HP:']
    if len(genes) <= 1:
        # print('Warning, only one gene in the community')
        return []
    tissues = []
    # find the tissues that gene is specific for
    tissue_specific_genes_count = 0
    for g in genes:
        temp_tissues = list(tissue_df[tissue_df['Description'] == g]['tissue'])
        if len(temp_tissues) > 0:
            tissue_specific_genes_count += 1
        tissues += temp_tissues

    if int(tissue_specific_genes_count * .5) <= 1:
        # print('Only one tissues specific gene')
        return []

    common_tissues = []
    for t in list(set(tissues)):
        if tissues.count(t) >= int(tissue_specific_genes_count * .5):
            common_tissues.append(t)
    return common_tissues


for comm in beckett_communities:
    tis = get_tissues_specific_to_all_genes(tsg,comm)
    if len(tis) > 0:
        print(tis)

print('\n')

half_count = 0
com_index = 0
for comm in beckett_communities:
    tis = get_tissues_specific_to_half_genes(tsg,comm)
    if len(tis) > 0:
        print(com_index)
        half_count += 1
        print(tis)
    com_index += 1
print(half_count)




