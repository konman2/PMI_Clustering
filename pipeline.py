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
name = 'wiki-fields'
filename = f'Predictions/{name}/new/predicted_communities_{taus}.npy'
predictions = np.load(filename)
file = f'Graphs/{name}/network.pkl'



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

comp = np.load(f'./trans_{name}_{len(times)}.npy')
# outfile = open(f'./results/pmi_values_{name}', 'w')
results_file = f'./results/results_{name}-(0,1000,50).txt'
results_file_strict = f'./results/results_{name}_strict-(0,1000,50).txt'
results = open(results_file).readlines()
results_strict = open(results_file_strict).readlines()

opt_pmi_strict = float(results_strict[1][:-1])
micro_nmi_strict = float(results_strict[-3][:-1])
macro_nmi_strict = float(results_strict[-2][:-1])
opt_pmi = float(results[1][:-1])
micro_nmi = float(results[-3][:-1])
macro_nmi = float(results[-2][:-1])
print(opt_pmi,opt_pmi_strict)
pmis = []
nmis_micro = []
nmis_macro = []
m = 0
loc = 0
for c,p in enumerate(predictions):
    val = calc_pmi(comp,p)[0]
    print(c,val)
    pmis.append(val)
    nmis_micro.append(normalized_mutual_info_score(p,gt_micro))
    nmis_macro.append(normalized_mutual_info_score(p,gt_macro))
nmis_micro = np.array(nmis_micro)
nmis_macro=np.array(nmis_macro)
times = np.linspace(-3,3,50)
# pmis=np.array(pmis)

fig1,ax1 = plt.subplots()
ax1.set_xlabel('Markov time log scale')
ax1.set_ylabel('PMI')
ax1.set_title(f'{name} PMIs at Different Times')
ax1.plot(times,pmis,label='varying')
ax1.plot(times,np.array([opt_pmi for i in times]),label='regular')
ax1.plot(times,np.array([opt_pmi_strict for i in times]),label='strict')
ax1.legend()
fig1.savefig(f'./results/plot_{name}_pmi')

fig2,ax2 = plt.subplots()
ax2.set_xlabel('Markov time log scale')
ax2.set_ylabel('NMI')
ax2.set_title(f'{name} NMI at Different Times (Micro)')
ax2.plot(times,nmis_micro,label='varying')
ax2.plot(times,np.array([ micro_nmi for i in times]),label='regular')
ax2.plot(times,np.array([ micro_nmi_strict for i in times]),label='strict')
ax2.legend()
fig2.savefig(f'./results/plot_{name}_nmi_micro')

fig3,ax3 = plt.subplots()
ax3.set_xlabel('Markov time log scale')
ax3.set_ylabel('NMI')
ax3.set_title(f'{name} NMI at Different Times (Macro)')
ax3.plot(times,nmis_macro,label='varying')
ax3.plot(times,np.array([ macro_nmi for i in times]),label='regular')
ax3.plot(times,np.array([ macro_nmi_strict for i in times]),label='strict')
ax3.legend()
fig3.savefig(f'./results/plot_{name}_nmi_macro')

#     print(c,val, normalized_mutual_info_score(p,gt_micro),normalized_mutual_info_score(p,gt_macro),file=outfile)
#     if val > m:
#         loc = c
#         m = val
# print("MAX: ",loc,m,file=outfile)