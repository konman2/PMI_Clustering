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
    #clust_a = [[i] for i in range(Mnn.shape[0])]
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
            #print("Node:"v)
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



def find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi,PMI):
    Graphs = []
    #PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
    clusters = np.arange(P.shape[0])
    hierch_clusters = []
    increased = True 
    next_clusters = np.zeros(P.shape[0])
    valid_locs = np.arange(P.shape[0])
    s = time.time()
    increased,n_Mnn,n_Mcn,n_Mnc,n_clusters,PMI = onestep(PMI,clusters,0.0001,Mcc,Mcn,Mnn,Mnc,Pi)
    e =  time.time()
    #print(e-s)
    #print(n_clusters)
    uniq_vals,next_clusters = np.unique(n_clusters,return_inverse=True)
    next_Mnn = n_Mnn[uniq_vals][:,uniq_vals]
    #print(next_clusters)
    #print(PMI)
    #exit()
    hierch_clusters.append(next_clusters)
    count = 1
    #print(count)
    
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
    for v in range(Mnn.shape[0]):
        c = v
        for p in hierch_clusters:
            c = p[c]
        final_cluster[v] = c
    return (final_cluster,PMI)
def run(flag='p',file='Graphs/airport_ww/network.pkl',name='airport',times='micro',precomp=None):
    #G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
    #G = nx.karate_club_graph()
    G = nx.read_gpickle(file)
    folder = 'pmi'
    if flag == 'ac':
        folder = 'ac'
    elif flag == 'ma':
        folder= 'mac'
    elif flag == 'lp':
        folder = 'lmepmi'
    print(G.number_of_nodes())
    #print(f.degree_vector(G).shape)
    #print(G.nodes)


    d = f.degree_vector(G)
    vert = np.arange(len(d))
    #print(d)
    # D= np.diag(1 / f.degree_vector(G))
    # A = f.adjacency_matrix(G)
    #P = np.asarray(np.matmul(D, A))
    #P = f.pagerank_transition_matrix(G)
    P = f.standard_random_walk_transition_matrix(G)
    Pi = d/(np.sum(d))
    #print(Pi)
 
    a = np.arange(10)+1
    y0 =[]
    y1= []
    y2 = []
    if times == 'macro' and (flag == 'lp' or flag == 'ma'):
        times == np.linspace(1,50,50)
    else:
        times = 10**np.linspace(-2,2,50)
    
    
    P_orig = np.copy(P)
    #P_orig = P_orig.astype('float128')
    #print(P)
    predictions = []
    time_taken = []
    num_clusters = []
    P_cum = np.zeros(P_orig.shape)
    P = None
    comp = []
    final_pmis = []
    for iters,t in enumerate(times):
        print("Run:",iters+1,"Markov Time:",t)
        s = time.time()
        #P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
        if precomp == None:
            print(P_orig.dtype)
            e_mat = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
           
            #test = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0])))
            print(np.min(e_mat[np.nonzero(e_mat)]))
            print(e_mat.dtype)
            print(np.argwhere(e_mat==0))
            comp.append(e_mat)
        else:
            e_mat = precomp[iters]
        if flag == 'lp' or flag == 'ma':
            #P_cum+=sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
            P_cum += e_mat
            P = P_cum/(iters+1)
        else:
            #P_cum = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
            P_cum = e_mat
            P = P_cum
        if flag == 'ac' or flag == 'ma':
            PMI = np.sum((Pi*P.diagonal()))-np.sum((Pi**2))
        else:
            PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
        #print("INITIIAL:",PMI)

        next_PMI = -1
        Mnn = None
        if flag == 'ac' or flag == 'ma':
            Mnn = ((P-Pi).T*Pi).T
        else:
            Mnn = np.log(P)-np.log(Pi)
        Mcc = np.copy(Mnn)
        Mcn = np.copy(Mnn)
        Mnc = np.copy(Mnn)
    
        final_cluster,pmi = find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi,PMI)
        final_pmis.append(pmi)
        print(pmi)
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
        #np.save('Predictions/airport/pmi/predicted_communities_{}'.format((iters+1)),np.array(predictions).T)
    #print(mapped_clusters)
    np.save('Predictions/{}/{}/predicted_communities_{}'.format(name,folder,len(times)),np.array(predictions).T)
    np.save('Predictions/{}/{}/times'.format(name,folder),times)
    np.save('Predictions/{}/{}/time_taken'.format(name,folder),np.array(time_taken))
    np.save('Predictions/{}/{}/num_clusters'.format(name,folder),np.array(num_clusters))
    np.save('Predictions/{}/{}/values'.format(name,folder),np.array(final_pmis))
    return comp

print('PMI')
# name = 'entsoe'
# file = 'Graphs/entsoe/network.pkl'
# name = 'airport'
# file = 'Graphs/airport_ww/network.pkl'
comp = list(np.load('computed_airport.npy'))
run('p',file=file,name=name,times='micro',precomp=None)
# np.save('computed_{}'.format(name),comp)
# print()
# print('LMEPMI')
# run('lp',file=file,name=name,times='micro',precomp=comp)
# print()
# print('AC')
#run('ac',file=file,name=name,times='micro',precomp=comp)
# print()
# print('MAC')
# run('ma',file=file,name=name,times='micro',precomp = comp)
