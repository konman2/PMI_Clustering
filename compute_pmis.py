import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time
import copy
import sys
import os
from make_custom import *
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score,adjusted_mutual_info_score
import pandas as pd

def calc_pmi(comp,clust,locs=None):
    if str(type(locs)) == '<class \'NoneType\'>':
        locs = np.ones(len(clust), dtype=bool)

    uniq_curr,cc = np.unique(clust,return_inverse=True)
    total = 0
    clust_scores = []
    for c in uniq_curr:
        Mnn = comp[0]
        clust_locs = (clust==c) & locs
        vals = np.sum(np.sum(comp[:,clust_locs][:,:,clust_locs],axis=1),axis=1)
        t = np.argmax(vals)
        best = vals[t]
        best_time = t
   
        total+=best
        clust_scores.append((c,best,best_time))
    return (total,clust_scores)
taus=50
name = 'airport_ww'
filename = f'Predictions/{name}/new/predicted_communities_{taus}.npy'
predictions = np.load(filename)
file = 'Graphs/airport_ww/network.pkl'



G = nx.read_gpickle(file)
gt_micro = np.load(f'Graphs/{name}/micro_comms.npy')
gt_macro = np.load(f'Graphs/{name}/macro_comms.npy')
# if len(subgraph_locs)> 0 :
#     G = G.subgraph(subgraph_locs+1)
#     G = nx.relabel.convert_node_labels_to_integers(G)
#     gt_micro = gt_micro[subgraph_locs]
#     gt_macro = gt_macro[subgraph_locs]



sigma = 1e-13
d = f.degree_vector(G)
vert = np.arange(len(d))

P = []
if nx.is_connected(G):
    P = f.standard_random_walk_transition_matrix(G)
else:
    P = f.pagerank_transition_matrix(G,mu=0.1)
Pi = d/(np.sum(d))

a = np.arange(10)+1
y0 =[]
y1= []
y2 = []

times = 10**np.linspace(-3,3,50)

#times = np.linspace(10**(-3),10**3,2000)
#times = np.linspace(10**(-3),10,200)
P_orig = np.copy(P)
#P_orig = P_orig.astype('float128')
#print(P)

time_taken = []
num_clusters = []
P_cum = np.zeros(P_orig.shape)
P = None
comp = []
final_pmis = []
runtimes = []
s= time.time()
# for iters,t in enumerate(times):
#     #print(t)
#     P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
#     M = np.log(P+sigma)-np.log(Pi+sigma)
#     comp.append(M)
#     print(t)
# comp = np.array(comp)
# np.save(f'./trans_airport_ww_{len(times)}',comp)
name = 'entsoe'
comp = np.load(f'./trans_{name}_{len(times)}.npy')
outfile = open(f'./results/pmi_values_{name}', 'w')
print(o,len(_),file=outfile)
m = 0
loc = 0
for c,p in enumerate(predictions):
    val = calc_pmi(comp,p)[0]
    print(c,val, normalized_mutual_info_score(p,gt_micro),normalized_mutual_info_score(p,gt_macro),file=outfile)
    if val > m:
        loc = c
        m = val
print("MAX: "c,val,file=outfile)