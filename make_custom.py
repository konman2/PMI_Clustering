import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time
import copy
import sys
import os
import itertools

def make_clique_tree(n,num_cliques):
    G = nx.Graph()
    for i in range(num_cliques):
        new_nodes = [j for j in range(i*n+1,i*n+n+1)]
        G.add_nodes_from(new_nodes)
        G.add_edges_from(list(itertools.combinations(new_nodes,2)))
        
    for i in range(num_cliques):
        if num_cliques<=1:
            break
        # if i == 0 and num_cliques > 1:
        #     G.remove_edge(1,n)
        #     G.add_edge(n,n+1)
        elif i == num_cliques-1:
            G.remove_edge(i*n+1,i*n+n)
        else:
            G.remove_edge(i*n+1,i*n+n)
            G.add_edge(i*n+n,i*n+n+1)
    G.add_edge(1,n*num_cliques)
    clust = np.array([(node-1)//n for node in G])
    return (G,clust,clust)
def asymmetric(arr,num_cliques):
    G = nx.Graph()
    curr_count = 1
    for i,n in enumerate(arr):
        new_nodes = [j for j in range(curr_count,curr_count+n)]
        print(new_nodes)
        
        #new_nodes = [j for j in range(i*n+1,i*n+n+1)]
        G.add_nodes_from(new_nodes)
        G.add_edges_from(list(itertools.combinations(new_nodes,2)))
        curr_count+=n
    print(G.nodes)
    cc = 1
    for i in range(num_cliques):
        n = arr[i]
        if num_cliques<=1:
            break
        # if i == 0 and num_cliques > 1:
        #     G.remove_edge(1,n)
        #     G.add_edge(n,n+1)
        elif i == num_cliques-1:
            G.remove_edge(cc,cc+n-1)
        else:
            G.remove_edge(cc,cc+n-1)
            G.add_edge(cc+n-1,cc+n)
        cc += n
    # G.add_edge(1,cc-1)
    # print(cc-1)
    nums = [[i for c in range(val)] for i,val in enumerate(arr)]
    clust = []
    for a in nums:
        clust+=a
    
    clust = np.array(clust)
    print(clust)
    return (G,clust,clust)
def hiearch(n,num_cliques):
    G = nx.Graph()
    cliques = []
    for i in range(num_cliques*(n)):
        new_nodes = [j for j in range(i*n+1,i*n+n+1)]
        G.add_nodes_from(new_nodes)
        G.add_edges_from(list(itertools.combinations(new_nodes,2)))
        cliques.append( new_nodes)
    print(cliques)
    last = 0
    for i in range(num_cliques):
        for j in range(last,last+n):
            for k in range(j+1,last+n):
                G.add_edges_from([(cliques[j][v],cliques[k][v]) for v in range(len(cliques[j]))])
        last+=n
    for i in range(num_cliques):
        if i != num_cliques-1:
            G.add_edge(cliques[(i+1)*n-1][-1],cliques[(i+1)*n][0])
    
            #print(last)
    print(G.nodes)
    print(G.edges)
    micro_clust = np.array([(node-1)//n for node in G])
    macro_clust = np.array([(node-1)//(n*n) for node in G])
    return (G,micro_clust,macro_clust)
def hiearch_weighted(n,num_cliques,w=1):
    G = nx.Graph()
    cliques = []
    for i in range(num_cliques*(n)):
        new_nodes = [j for j in range(i*n+1,i*n+n+1)]
        G.add_nodes_from(new_nodes)
        G.add_edges_from(list(itertools.combinations(new_nodes,2)))
        cliques.append( new_nodes)
    print(cliques)
    last = 0
    for i in range(num_cliques):
        for j in range(last,last+n):
            for k in range(j+1,last+n):
                new_edges = []
                for v in range(len(cliques[j])):
                    if np.random.binomial(1,w) == 1:
                        new_edges.append((cliques[j][v],cliques[k][v]))
                G.add_edges_from(new_edges)
                #G.add_edges_from([(cliques[j][v],cliques[k][v]) for v in range(len(cliques[j]))])
        last+=n
    for i in range(num_cliques):
        if i != num_cliques-1:
            G.add_edge(cliques[(i+1)*n-1][-1],cliques[(i+1)*n][0])
    
            #print(last)
    print(G.nodes)
    print(G.edges)
    micro_clust = np.array([(node-1)//n for node in G])
    macro_clust = np.array([(node-1)//(n*n) for node in G])
    return (G,micro_clust,macro_clust)

#G,micro_comms,macro_comms = hiearch(6,5)
# n = 30
#arr = [30,10,20,50,5]
#G,micro_comms,macro_comms = make_clique_tree(30,3)
#G,micro_comms,macro_comms = hiearch_weighted(6,5,0.75)
#G,micro_comms,macro_comms = hiearch_weighted(10,5,0.05)

# G,micro_comms,macro_comms = make_clique_tree(n,10)
#G,clust = asymmetric(arr,len(arr))
# clust = np.array([(node-1)//n for node in G])
# micro_comms = clust
#print(micro_comms)
# np.save('Graphs/custom/micro_comms.npy',micro_comms)
# np.save('Graphs/custom/macro_comms.npy',macro_comms)
# nx.write_gpickle(G,'Graphs/custom/network.pkl')