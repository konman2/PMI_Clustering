import networkx as nx
G = nx.read_gpickle('Graphs/snp500ll/alpha=0.55/network.pkl')
print(G.nodes)
print(nx.algorithms.components.number_connected_components(G))
print(nx.bipartite.is_bipartite(G))