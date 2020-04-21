import networkx as nx
from Becketts_Community_Analysis_Results.becketts_community_analysis import get_beckett_communities, load_graphs, webweb_plot


# load the MyGene2 Networks
G, Gn = load_graphs()

# find communities using Greedy
c = list(nx.algorithms.community.modularity_max.greedy_modularity_communities(G))
# convert from 2D frozenset to 2D list
greedy_coms = [[y for y in x] for x in c]
# plot greedy
webweb_plot(Gn, greedy_coms)

# load the Beckett communities
beckett_coms = get_beckett_communities()
# plot Beckett
webweb_plot(Gn, beckett_coms)
