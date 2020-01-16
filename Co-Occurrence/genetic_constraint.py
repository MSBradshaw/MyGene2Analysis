import pandas as pd
import numpy as np
import pickle

# read in the tabular data
data = pd.read_csv('Data/my_gene_2_variantes_by_family_tables.csv')

# get the number of genes and number of hpos
genes = list(set(data.iloc[:, 1]))

lof = pd.read_csv('gnomad.v2.1.1.lof_metrics.by_gene.txt', sep='\t')
# example
lof[lof['gene'] == 'SS18L1'].loc[:, 'oe_lof_lower']

pLI = []
oe_lof = []
oe_lof_lower = []
oe_lof_upper = []
for g in genes:
    p = lof[lof['gene'] == g].loc[:, 'pLI'].values
    ol = lof[lof['gene'] == g].loc[:, 'oe_lof'].values
    oll = lof[lof['gene'] == g].loc[:, 'oe_lof_lower'].values
    olu = lof[lof['gene'] == g].loc[:, 'oe_lof_upper'].values
    if len(p) > 0:
        pLI.append(p[0])
    else:
        pLI.append(p)
    if len(ol) > 0:
        oe_lof.append(ol[0])
    else:
        oe_lof.append(ol)
    if len(oll) > 0:
        oe_lof_lower.append(oll[0])
    else:
        oe_lof_lower.append(oll)
    if len(olu) > 0:
        oe_lof_upper.append(olu[0])
    else:
        oe_lof_upper.append(olu)

df = pd.DataFrame(
    {'gene': genes, 'pLI': pLI, 'oe_lof': oe_lof, 'oe_lof_lower': oe_lof_lower, 'oe_lof_upper': oe_lof_upper})
df.to_csv('Data/gene-constraints.csv')

# cross reference the 26 candidates from the phi analysis with lof information
phi_sums_df = pickle.load(open("Data/phi_sums_df.pickle", "rb"))

phi_sums_small = phi_sums_df[phi_sums_df['hpo_sum'] > 2]
phi_sums_small = phi_sums_small[phi_sums_small['phi'] > 0.3]

phi_sums_small = phi_sums_small.sort_values(by=['phi'],ascending=False)

df.index = df['gene']
phi_sums_small.index = phi_sums_small['Gene']
df.intersection(phi_sums_small)

sub_set_lof = df.loc[df.index.intersection(phi_sums_small.index),:]
sub_set_lof.to_csv('Data/gene-constraints-subset.csv')