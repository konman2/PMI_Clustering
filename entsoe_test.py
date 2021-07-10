import networkx as nx
G = nx.read_gpickle('Graphs/entsoe/network.pkl')
print(nx.algorithms.components.number_connected_components(G))
print(nx.bipartite.is_bipartite(G))