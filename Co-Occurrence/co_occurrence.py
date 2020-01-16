import pandas as pd
import numpy as np
import pickle
import os.path
from os import path
import matplotlib.pyplot as plt
import seaborn as sns
from math import sqrt
import re


def calc_phi(tp, fp, tn, fn):
    # https://stackoverflow.com/a/56875660/992687
    x = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    return ((tp * tn) - (fp * fn)) / sqrt(x)


co = None
if path.exists("Data/co_occurrence_pandas.pickle"):
    print('Loading Pre-made co-occurrence matrix')
    co = pickle.load(open("Data/co_occurrence_pandas.pickle", "rb"))
else:
    print('Generating co-occurrence matrix from scratch')
    data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

    # get the number of genes and number of hpos
    genes = list(set(data.iloc[:, 1]))

    # the number of hpos and hpo names do not match. There are more names than hpo IDs
    hpos = list(set(data.iloc[:, 7]))

    # create an empty co-occurance matrix with genes as rows and hpos as columns
    co = pd.DataFrame(np.zeros((len(genes), len(hpos))))
    co.columns = hpos
    co.index = genes

    for i in range(data.shape[0]):
        row = data.iloc[i, :]
        co.loc[row[1], row[7]] += 1
    print('Pickling co-occurrence matrix for future use')
    pickle.dump(co, open("Data/co_occurrence_pandas.pickle", "wb"))

# remove row and columns with a sum of 20
rows_drop = []
row_sums = []
for i in co.index:
    row_sums.append((i, sum(co.loc[i, :])))
    if sum(co.loc[i, :]) < 100:
        rows_drop.append(i)

# remove row and columns with a sum of 20
col_drop = []
col_sums = []
for i in co.columns:
    col_sums.append((i, sum(co.loc[:, i])))
    if sum(co.loc[:, i]) < 50:
        col_drop.append(i)

row_sums.sort(key=lambda x: x[1])
row_sums.reverse()
col_sums.sort(key=lambda x: x[1])
col_sums.reverse()

peaks = []
for i in co.index:
    for j in co.columns:
        peaks.append((i + '-' + j, co.loc[i, j]))
peaks.sort(key=lambda x: x[1])
peaks.reverse()

for i in range(20):
    print(peaks[i])

print_co = co.copy()
# add 1 to the columns with the highest sums so we know which are which
# for i in range(6):
#     print(col_sums[i][0])
#     print_co.loc[:, col_sums[i][0]] = print_co.loc[:, col_sums[i][0]] + 5
#
# for i in range(10):
#     print(row_sums[i][0])
#     print_co.loc[row_sums[i][0], :] = print_co.loc[row_sums[i][0], :] + 5

co_small = print_co.drop(col_drop, 1)
co_small = co_small.drop(rows_drop, 0)
# plt.imshow(co_small.to_numpy())

# sns.set(rc={'figure.figsize': (18, 10)})
# ax = sns.heatmap(co_small, xticklabels=True, yticklabels=True)
# ax.figure.savefig('heatmap.png')

phi = None
if path.exists("Data/phi_coefficinets_pandas.pickle"):
    print('Loading Pre-made phi_coefficinets matrix')
    phi = pickle.load(open("Data/phi_coefficinets_pandas.pickle", "rb"))
else:
    print('Creating phi_coefficinets matrix from scratch')
    # make a blank data frame for recording the phi-coefficient
    phi = pd.DataFrame(np.zeros(co.shape))
    phi.columns = co.columns
    phi.index = co.index
    row_sums2 = {}
    for gene in co.index:
        row_sums2[gene] = False
    sum_of_co = sum(sum(co.to_numpy()))
    for hpo in co.columns:
        col_sum = sum(co.loc[:, hpo])
        for gene in co.index:
            # if there is only one mention of the gene, ignore it
            if col_sum == 1:
                phi.loc[gene, hpo] = -1
                continue
            true_pos = co.loc[gene, hpo]
            false_neg = 0
            # check if the row sum has already been calculated
            if row_sums2[gene] != False:
                false_neg = row_sums2[gene] - true_pos
            else:
                rs = sum(co.loc[gene, :])
                false_neg = rs - true_pos
                row_sums2[gene] = rs
            false_pos = col_sum - true_pos
            true_neg = sum_of_co - sum(co.loc[gene, :]) - sum(co.loc[:, hpo]) + true_pos
            p = calc_phi(true_pos, false_pos, true_neg, false_neg)
            phi.loc[gene, hpo] = p
    pickle.dump(phi, open("Data/phi_coefficinets_pandas.pickle", "wb"))

# what interesting information is there in the phi matrix
phi_list = []
for gene in phi.index:
    for hpo in phi.columns:
        phi_list.append((gene + '-' + hpo, phi.loc[gene, hpo]))
phi_list.sort(key=lambda x: x[1])
phi_list.reverse()

phi_sums_df = None
if path.exists("Data/phi_sums_df.pickle"):
    print('Loading Pre-made phi_sums_df')
    phi_sums_df = pickle.load(open("Data/phi_sums_df.pickle", "rb"))
else:
    print('creating phi_sums_df from scratch')
    phi_sums_df = pd.DataFrame({'Gene': [], 'HPO': [], 'phi': [], 'gene_sum': [], 'hpo_sum': []})
    gene_sums = {}
    hpo_sums = {}

    for i in co.index:
        gene_sums[i] = False

    for i in co.columns:
        hpo_sums[i] = False

    for i in range(len(phi_list)):
        g = re.sub(r'-HP.*', '', phi_list[i][0])
        h = re.sub(r'.*-', '', phi_list[i][0])
        if g == 'gene' or g == 'gene-phenotype_hpo_id' or 'gene' in g or 'gene' in h or '-pheno' in h or '-pheno' in g:
            continue
        gs = -1
        hs = -1
        if hpo_sums[h] != False:
            hs = hpo_sums[h]
        else:
            hs = sum(co.loc[:, h])
            hpo_sums[h] = hs
        if gene_sums[g] != False:
            gs = gene_sums[g]
        else:
            gs = sum(co.loc[g, :])
            gene_sums[g] = gs

        p = phi_list[i][1]
        phi_sums_df = phi_sums_df.append({'Gene': g, 'HPO': h, 'phi': p, 'gene_sum': gs, 'hpo_sum': hs}, ignore_index=True)
    pickle.dump(phi_sums_df, open("Data/phi_sums_df.pickle", "wb"))


# sort the df by phi

phi_sums_small = phi_sums_df[phi_sums_df['hpo_sum'] > 2]
phi_sums_small = phi_sums_small[phi_sums_small['phi'] > 0.3]

phi_sums_small = phi_sums_small.sort_values(by=['phi'],ascending=False)