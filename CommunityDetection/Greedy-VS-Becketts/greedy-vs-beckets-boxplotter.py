import pickle
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from os import path

plt.margins(x=0.1)


"""
Return of list of lists of communities
"""


def get_greedy_communities():
    greedy_pickle_path = 'CachedFiles/greedy_communities.pickle'
    if path.exists(greedy_pickle_path):
        print('Loading Greedy Communities')
        return pickle.load(open(greedy_pickle_path, 'rb'))
    else:
        print('Generating Greedy Communities')
        G, Gn = load_graphs()
        c = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(G))
        # convert from 2D frozenset to 2D list
        return [[y for y in x] for x in c]


# use function from other file to load the Becketts communities
communities = get_beckett_communities()

becketts_lengths = []
for c in communities:
    becketts_lengths.append(len(c))

g_coms = get_greedy_communities()

greedy_lengths = []
for c in g_coms:
    greedy_lengths.append(len(c))

com_lengths = np.array(greedy_lengths + becketts_lengths)
com_lengths_log = np.log(com_lengths)

# Boxplot raw data
greedy_name = 'Greedy\n N=' + str(len(greedy_lengths))
becketts_name = 'LPAwb+\n N=' + str(len(becketts_lengths))
data = pd.DataFrame({'Size of Communities': com_lengths,
                     'Algorithm': [greedy_name] * len(greedy_lengths) + [becketts_name] * len(becketts_lengths)})

ax1 = sns.boxplot(x="Algorithm", y="Size of Communities", data=data)
ax1.set_xlabel('')
ax1.figure.savefig("CommunityDetection/Greedy-VS-Becketts/greedy-v-beckett-box.png")
plt.clf()

# Violin Plot with Logged Data
greedy_name = 'Greedy\n N=' + str(len(greedy_lengths))
becketts_name = 'LPAwb+\n N=' + str(len(becketts_lengths))
data = pd.DataFrame({'Log(Size of Communities)': com_lengths_log,
                     'Algorithm': [greedy_name] * len(greedy_lengths) + [becketts_name] * len(becketts_lengths)})

ax2 = sns.violinplot(x="Algorithm", y="Log(Size of Communities)", data=data, split=True,
                     scale="count", inner="stick")
ax2.set_xlabel('')
ax2.figure.savefig("CommunityDetection/Greedy-VS-Becketts/greedy-v-beckett-log-violin-log.png")
plt.clf()

# Violin Plot with raw data minus greedy super community
greedy_name = 'Greedy\n N=' + str(len(greedy_lengths) - 1)
becketts_name = 'LPAwb+\n N=' + str(len(becketts_lengths))
data = pd.DataFrame({'Size of Communities': com_lengths[1:],
                     'Algorithm': [greedy_name] * (len(greedy_lengths) - 1) + [becketts_name] * len(becketts_lengths)})

ax3 = sns.violinplot(x="Algorithm", y="Size of Communities", data=data, split=True,
                     scale="count", inner="stick")
ax3.set_xlabel('')
ax3.figure.savefig("CommunityDetection/Greedy-VS-Becketts/greedy-v-beckett-outlier-excluded-violin.png")
plt.clf()
