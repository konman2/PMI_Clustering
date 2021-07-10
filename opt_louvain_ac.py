import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time
import copy
import sys

min_inc = 0.0001
def onestep(PMI,clusters,min_inc,Mcc,Mcn,Mnn,Mnc,Pi):
    tol =1e-6
    increased = False
    next_PMI = -1
    count = 0
  
    while (next_PMI-PMI > min_inc or next_PMI == -1) and count<1000:
        if next_PMI == -1:
            next_PMI = PMI
        else:
            PMI = next_PMI
       
        s = time.time()
        
        for v in range(Mnn.shape[0]):
            best=v
            best_inc = 0
            cv = clusters[v]
            decr = -Mcc[cv,cv] 
            Pcvcv_new = (Mcc[cv][cv]-Mnc[v][cv]-Mcn[cv][v]+Mnn[v][v])

            decr += Pcvcv_new
          
            inc_pmi = Mcn[:,v]+Mnc[v]+Mnn[v][v]+decr
            inc_pmi[cv] = 0
            best = np.argmax(inc_pmi)
            best_inc = inc_pmi[best]
            if best_inc <= tol or best == clusters[v]:
                    continue
            else:
                increased = True
            c = best
            cv = clusters[v]
            #updates
            Mcc_new = ((Mcc[c][c]+Mcn[c][v]) +(Mnc[v][c]+Mnn[v][v]))
            Pcvcv_new = (Mcc[cv][cv]-Mnc[v][cv]-Mcn[cv][v]+Mnn[v][v])
            Pccv_new = Mcc[c][cv]+Mnc[v][cv]-Mnn[v][v]-Mcn[c][v]
            Pcvc_new = Mcc[cv][c]-Mnc[v][c]+Mcn[cv][v]-Mnn[v][v]
            Mcc[c]+=Mnc[v]
            Mcc[cv]-=Mnc[v]
            Mcc[:,c]+=Mcn[:,v]
            Mcc[:,cv]-=Mcn[:,v]
            Mcn[c]+=Mnn[v]
            Mcn[cv]-=Mnn[v]
            Mnc[:,c]+=Mnn[:,v]
            Mnc[:,cv]-=Mnn[:,v]   
            Mcc[c][c] = Mcc_new
            Mcc[cv][cv]= Pcvcv_new
            Mcc[c][cv]=Pccv_new
            Mcc[cv][c]=Pcvc_new
            clusters[v] = c
            next_PMI += best_inc
         
        e = time.time()
        count+=1
        print()
        print("Run:",count,"Time:",e-s,"PMI:",next_PMI)
        print()
    return (increased,Mcc,Mcn,Mnc,clusters,next_PMI)



def find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi):
    Graphs = []
    PMI = np.sum(Pi*P.diagonal())-np.sum(Pi**2)
    clusters = np.arange(P.shape[0])
    hierch_clusters = []
    increased = True 
    next_clusters = np.zeros(P.shape[0])
    valid_locs = np.arange(P.shape[0])
    s = time.time()
    increased,n_Mnn,n_Mcn,n_Mnc,n_clusters,PMI = onestep(PMI,clusters,0.0001,Mcc,Mcn,Mnn,Mnc,Pi)
    
    e =  time.time()

    uniq_vals,next_clusters = np.unique(n_clusters,return_inverse=True)
    next_Mnn = n_Mnn[uniq_vals][:,uniq_vals]
    #exit()
    hierch_clusters.append(next_clusters)
    count = 1
    while increased:
        increased = False
        Mnn_new = next_Mnn
        Mcn_new = np.copy(Mnn_new)
        Mnc_new = np.copy(Mnn_new)
        Mcc_new = np.copy(Mnn_new)
        new_clusters = np.arange(Mnn_new.shape[0])
        increased,next_Mnn,next_Mcn,next_Mnc,next_clusters,PMI = onestep(PMI,new_clusters,0.0001,Mcc_new,Mcn_new,Mnn_new,Mnc_new,Pi)
        uniq_vals,next_clusters = np.unique(next_clusters,return_inverse=True)
        next_Mnn = next_Mnn[uniq_vals][:,uniq_vals]
        hierch_clusters.append(next_clusters)
        count+=1
        print("Number of runs:", count)
        
    final_cluster = np.zeros(P.shape[0])
    #print(hierch_clusters)
    for v in range(Mnn.shape[0]):
        c = v
        for p in hierch_clusters:
            c = p[c]
        final_cluster[v] = c
    print(PMI)
    return final_cluster

#G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
#G_1 = nx.karate_club_graph()
G = nx.read_gpickle('Graphs/airport_ww/network.pkl')
print(G.number_of_nodes())
print(f.degree_vector(G).shape)
print(G.nodes)


d = f.degree_vector(G)
vert = np.arange(len(d))
print(d)
D= np.diag(1 / f.degree_vector(G))
A = f.adjacency_matrix(G)
P = np.asarray(np.matmul(D, A))

Pi = d/(np.sum(d))
print(Pi)
a = np.arange(10)+1
y0 =[]
y1= []
y2 = []
times = 10**np.linspace(-2,-2,100)

#times = [1]
#P = P.todense()
P_orig = np.copy(P)
print(P)
predictions = []
time_taken = []
num_clusters = []
for iters,t in enumerate(times):
    print("Run:",iters+1,"Markov Time:",t)
    s = time.time()
    P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)

    PMI = np.sum((Pi*P.diagonal()))-np.sum((Pi**2))

    #print("INITIIAL:",PMI)

    next_PMI = -1
    Mnn = ((P-Pi).T*Pi).T
    Mcc = np.copy(Mnn)
    Mcn = np.copy(Mnn)
    Mnc = np.copy(Mnn)
   
    final_cluster = find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi)
    e = time.time()
    time_taken.append(e-s)
    num_clusters.append(len(np.unique(final_cluster)))
    map_clust = {}
    mapped_clusters = []
    count = 1
    for i in final_cluster:
        if i not in map_clust:
            map_clust[i] = count
            count+=1
        mapped_clusters.append(map_clust[i])
    predictions.append(mapped_clusters)
    
    #np.save('Predictions/ac/predicted_communities_{}'.format((iters+1)),np.array(predictions).T)
#print(mapped_clusters)
np.save('Predictions/ac/predicted_communities_{}'.format(len(times)),np.array(predictions).T)
np.save('Predictions/ac/times',times)
np.save('Predictions/ac/time_taken',np.array(time_taken))
np.save('Predictions/ac/num_clusters',np.array(num_clusters))
