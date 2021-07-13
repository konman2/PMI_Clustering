import networkx as nx
G = nx.read_gpickle('Graphs/wiki-fields/network.pkl')
print(nx.algorithms.components.number_connected_components(G))
print(nx.bipartite.is_bipartite(G))