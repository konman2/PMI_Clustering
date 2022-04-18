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
min_inc = 0.0001


def find_cluster(seed_node,min_inc,Mnn,debug=False):
    Mcn = Mnn[seed_node]
    Mnc = Mnn[:,seed_node]
    tol =1e-13
    increased = True
    count = 0
    clusters = np.ones(Mnn.shape[0])*-1
    clusters[seed_node] = 1
    if debug:
        print(Mcn)
        print(Mnc)
    while increased:
        increased= False
        s = time.time()
        inc_pmi = Mcn+Mnc
        inc_pmi[clusters == 1] = 0
        best = np.argmax(inc_pmi)
        best_inc = inc_pmi[best]
        #print(best,best_inc)
        if best_inc <= tol or clusters[best] == 1:
            continue
        else:
            increased = True
        clusters[best] = 1
        Mcn+= Mnn[best]
        Mnc+=Mnn[:,best]
    
    return clusters


def onestep(PMI,clusters,min_inc,Mcc,Mcn,Mnn,Mnc,Pi,sample=False,k=0):
    #clust_a = [[i] for i in range(Mnn.shape[0])]
    tol =1e-13
    increased = False
    next_PMI = -1
    count = 0
    while (next_PMI-PMI > tol or next_PMI == -1) and count<1000:
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
            #print(c,v)
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
        # print()
        # print("Run:",count,"Time:",e-s,"PMI:",next_PMI)
        # print()
    
    return (increased,Mcc,Mcn,Mnc,clusters,next_PMI)



def find_clusters(Mcc,Mnn,Mcn,Mnc,Pi,PMI):
    Graphs = []
    #PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
    clusters = np.arange(Mnn.shape[0])
    hierch_clusters = []
    increased = True 
    next_clusters = np.zeros(Mnn.shape[0])
    valid_locs = np.arange(Mnn.shape[0])
    s = time.time()
    increased,n_Mnn,n_Mcn,n_Mnc,n_clusters,PMI = onestep(PMI,clusters,0.0001,Mcc,Mcn,Mnn,Mnc,Pi)
    e =  time.time()
    uniq_vals,next_clusters = np.unique(n_clusters,return_inverse=True)
    next_Mnn = n_Mnn[uniq_vals][:,uniq_vals]
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
    final_cluster = np.zeros(Mnn.shape[0])
    for v in range(Mnn.shape[0]):
        c = v
        for p in hierch_clusters:
            c = p[c]
        final_cluster[v] = c
    #print(final_cluster)
    return (final_cluster,PMI)

def func(M,locs,nc,c,uniq_curr):
    Mnn = np.zeros((len(locs),len(locs)))
    for i in uniq_curr:
        for j in uniq_curr:
            Mnn[i][j] = np.sum(M[locs[i]][:,locs[j]])
    l = nc == c
    return np.sum(np.sum(Mnn[l][:,l]))

def find_best2(locs,comp,deg_vec,start=0):
    #print(comp.shape)
    vals = np.sum(np.sum(comp[:,locs][:,:,locs],axis=1),axis=1)
    t = np.argmax(vals)
    best = vals[t]
    # Mnn = comp[0]
    # time = 0
    # # best =  np.sum(Mnn[np.where(nc == c)][:,np.where(nc==c)])
    # best = np.sum(Mnn[locs][:,locs])
    # for i in range(1,len(comp)):
    #     Mnn = comp[i]
    #     val = np.sum(Mnn[locs][:,locs])
    #     if val > best:
    #         time = i
    #         best = val
    #     # best = max(best,val)
    return best/np.sum(locs)
    # if np.sum(locs) == 1:
    #     print(comp[t])
    #return (best/np.sum(locs),t)
def find_best(comp,nc,c,start=0):
    Mnn = comp[0]
    locs = nc == c
    time = 0
    # best =  np.sum(Mnn[np.where(nc == c)][:,np.where(nc==c)])
    best = np.sum(Mnn[locs][:,locs])
    for i in range(1,len(comp)):
        Mnn = comp[i]
        val = np.sum(Mnn[locs][:,locs])
        if val > best:
            time = i
            best = val
        # best = max(best,val)
    return (best,time,locs)


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
        #best = np.sum(Mnn[clust_locs][:,clust_locs])
        best_time = t
        # for t in range(1,len(comp)):
        #     Mnn = comp[t]
        #     val = np.sum(Mnn[clust_locs][:,clust_locs])
        #     if val > best:
        #         best = val
        #         best_time = t
        total+=best
        clust_scores.append((c,best,best_time))
    return (total,clust_scores)


def arlei_graph():
    sizes = [160, 80, 40, 20]
    p = 0.03
    # sizes = [160,20,160,40]
    probs = .5 * np.identity(len(sizes))

    probs[0,1] = p
    probs[1,0] = p

    probs[2,3] = 0.02
    probs[3,2] = 0.02
    probs[0,2] = p
    probs[2,0] = p

    probs[1,2] = p
    probs[2,1] = p

    probs[1,3] = p
    probs[3,1] = p

    probs[0,3] = p
    probs[3,0] = p

    sbm = nx.stochastic_block_model(sizes, probs, seed=0)
    print(sbm.nodes)
    print(type(sbm))
    count = 0
    gt = []
    for a in sizes:
        for i in range(a):
            gt.append(count)
        count+=1
    gt = np.array(gt)
    return sbm,gt,gt

def run(flag='p',file='Graphs/airport_ww/network.pkl',name='airport',times='micro',precomp=None,subgraph_locs=[]):
    #G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
    #G = nx.karate_club_graph()
    # gt_micro = []
    # gt_macro = []
    
    if name == 'airport':
        name = 'airport_ww'
    gt_micro = np.load(f'Graphs/{name}/micro_comms.npy')
    gt_macro = np.load(f'Graphs/{name}/macro_comms.npy')


    G = nx.read_gpickle(file)
    if len(subgraph_locs)> 0 :
        G = G.subgraph(subgraph_locs+1)
        G = nx.relabel.convert_node_labels_to_integers(G)
        gt_micro = gt_micro[subgraph_locs]
        gt_macro = gt_macro[subgraph_locs]

    # gt = np.load(f'Graphs/{name}/micro_comms.npy')
    # assert len(gt) == len(G.nodes)

    # G,gt_micro,gt_macro = arlei_graph()
    # name="a"
    #print(gt)
    # exit()
    # arr = [40,20,10]
  
    # G,gt,_ = asymmetric(arr,len(arr))
    # G = nx.read_gpickle(file)
    # gt = np.load(f'Graphs/{name}/micro_comms.npy')
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

    sigma = 1e-13
    d = f.degree_vector(G)
    vert = np.arange(len(d))
    #print(d)
    # D= np.diag(1 / f.degree_vector(G))
    # A = f.adjacency_matrix(G)
    #P = np.asarray(np.matmul(D, A))
    #P = f.pagerank_transition_matrix(G,mu=0.1)
    #print(np.nonzero(P==0))
    P = []
    if nx.is_connected(G):
        P = f.standard_random_walk_transition_matrix(G)
    else:
        P = f.pagerank_transition_matrix(G,mu=0.1)
    Pi = d/(np.sum(d))
    #pi = f.stationary_distribution(P)
    #Pi = np.diag(pi)
    #print(Pi)
    #np.save('Pis_{}'.format(name),Pi) 
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
    predictions = []
    time_taken = []
    num_clusters = []
    P_cum = np.zeros(P_orig.shape)
    P = None
    comp = []
    final_pmis = []
    runtimes = []
    s= time.time()
    for iters,t in enumerate(times):
        #print(t)
        P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
        M = np.log(P+sigma)-np.log(Pi+sigma)
        comp.append(M)
    comp = np.array(comp)
    # print(comp[29])
    # exit()
    np.save(f'./trans_{name}_{len(times)}',comp)
    e = time.time()
    runtimes.append(e-s)
    print(comp.shape)
    #exit()


    curr_clusters = np.arange(P.shape[0])
    uniq_clust = np.arange(P.shape[0])
    best_times = np.zeros(len(uniq_clust))
    
    best_pmi = np.copy(np.diag(comp[0]))
    print(comp[0])
    print(best_pmi)
    predictions = []
    # exit()
    sourceFile = open('out.txt', 'w')

    s = time.time()
    selected = np.random.choice(P.shape[0],P.shape[0]//2)
    for iters,t in enumerate(times):
        M = np.copy(comp[iters])
        for s in selected:
            # if iters == 40:
                
            #     l = find_cluster(s,1e-13,M,False)
            #     print(M[s]+M[:,s])
            #     print(np.sum(l == 1))
            # else:
            l = find_cluster(s,1e-13,M,False)
            predictions.append(l)
        
    

    off_limits = set()
    signal = np.zeros(P.shape[0])
    count = 0
    sourceFile.close()
    e = time.time()
    runtimes.append(e-s)
    s = time.time()
    # if not os.path.isdir(f'Predictions/{name}'):
    #     os.mkdir(f'Predictions/{name}')
    # if not os.path.isdir(f'Predictions/{name}/new'):
    #     os.mkdir(f'Predictions/{name}/new')
    # np.save('Predictions/{}/new/predicted_communities_{}'.format(name,len(times)),np.array(predictions))
    # # np.save('Predictions/{}/{}/times'.format(name,folder),times)
    # # np.save('Predictions/{}/{}/time_taken'.format(name,folder),np.array(time_taken))
    # # np.save('Predictions/{}/{}/num_clusters'.format(name,folder),np.array(num_clusters))
    # # np.save('Predictions/{}/{}/values'.format(name,folder),np.array(final_pmis))
    print("started")
    predictions.append(np.arange(P.shape[0]))
    predictions = np.array(predictions)
    #print(predictions)
    total = 0
    while(total<P.shape[0]):
        best_l = -1
        best_ind  = -1
        best = 0
        #proc = set()
        #valid_locs = np.array([])
        valid_locs = []
        og_pred = []
        # s = time.time()
        for i,clustering in enumerate(predictions):
            clusters = np.unique(clustering)
            #print(clustering)
            #print
            for c in clusters:
                if c >= P.shape[0] or c<0:
                    continue
                l = clustering == c
                #z = str(l)
                l = l.reshape(1,P.shape[0])
                #print(len(valid_locs))
                if len(valid_locs) == 0:
                    valid_locs = np.copy(l)
                    og_pred.append(i)
                elif np.max(np.sum(valid_locs == l,axis=1)) < P.shape[0]:
                    valid_locs = np.append(valid_locs,l,axis=0)
                    og_pred.append(i)
                

                # if z not in proc:
                #     valid_locs.append(l)
                #     proc.add(z)
                
                #print(clustering,c)
        # e = time.time()
        # print(e-s)
        # print(len(valid_locs))
        # exit()
        # s = time.time()
        # for loc in valid_locs:
        #     val= find_best2(loc,comp,start=0)
        #     #val/=np.sum(locs)
        #     #print(val,i,t,c)
        #     #val/=(len(locs)**2)
        #     if val > best:
        #         best = val
        #         best_l = loc
                #print(best_c)
                #best_ind = i
        #print(valid_locs.shape)
        print("#locs",len(valid_locs))
        vals = np.apply_along_axis(find_best2,1,valid_locs,comp,d)
        # for c in range(len(valid_locs)):
        #     val,t = find_best2(valid_locs[c],comp,d)
        #     if np.sum(valid_locs[c]) == 1:
        #         print(comp[t])
        #     print(val,t,np.sum(valid_locs[c]))
        # exit()
        #print(vals)
        # exit()
        best_i = np.argmax(vals)
        print(best_i,vals.shape)
        best = vals[best_i]
        # assert np.sum(best_l == valid_locs[best_i]) == P.shape[0]
        best_l = valid_locs[best_i]
        # e = time.time()
        # print(e-s)
       
        next_clust = P.shape[0]+count
        # print(best_ind)
        #print(predictions[best_ind],best_c)
        #locs = (predictions[best_ind] == best_c)
        #print(np.where(locs),best,best_c,best_ind)
        
        predictions[:,best_l] = next_clust
        # for i in range(len(predictions)):
        #     predictions[i][best_l] != next_clust
        # off_limits.add(next_clust)
        count+=1
        signal[best_l] = 1
        total = np.sum(signal)
        print(best,np.sum(best_l),og_pred[best_i])
        print(total)
        #print(np.sum(best_l))
        # e = time.time()
        # print(e-s)
        # print()
        #print(np.sum(signal))
        #print(best,best_l)
        #print((np.sum(l),i) for i,l in enumerate(locs))
        #print(predictions[0])
    _,o = np.unique(predictions[0],return_inverse=True)
    e = time.time()
    runtimes.append(e-s)
    to_file=True
    if to_file:
        outfile = open(f'./results/results_{name}-({int(times[0])},{int(times[-1])},{len(times)})-seeds.txt', 'w')
        print(o,len(_),file=outfile)
        # z = np.copy(o)
        # z[z==2] = 1
        #print(z)
        #print(calc_pmi(comp,z)[0])
        print(calc_pmi(comp,o)[0],file=outfile)
        #print(gt)
        print(calc_pmi(comp,gt_micro)[0],file=outfile)
        print("",file=outfile)
        print(calc_pmi(comp,gt_macro)[0],file=outfile)
        print("",file=outfile)
        print(normalized_mutual_info_score(o,gt_micro),file=outfile)
        print(normalized_mutual_info_score(o,gt_macro),file=outfile)
        print(runtimes,file=outfile)
        outfile.close()
    else:
        print(o,len(_))
        # z = np.copy(o)
        # z[z==2] = 1
        #print(z)
        #print(calc_pmi(comp,z)[0])
        print(calc_pmi(comp,o)[0])
        #print(gt)
        print(calc_pmi(comp,gt_micro)[0])
        print("")
        print(calc_pmi(comp,gt_macro)[0])
        print("")
        print(normalized_mutual_info_score(o,gt_micro))
        print(normalized_mutual_info_score(o,gt_macro))
        print(runtimes)

  



            


        
    #     print(t)
    #     # print(nc,"fakjf")
    #     print(curr_clusters)
    #     # print(uniq_clust)
    #     test_clust = np.copy(curr_clusters)
    #     for c,i in enumerate(uniq_clust):
    #         test_clust[curr_clusters == c] = nc[i]
    #     print(test_clust)
    #     # print(test_clust == nc)
    #     #flattened_vals = 
    #     #clust = np.copy(curr_clusters)
    #     temp = np.copy(curr_clusters)
    #     for c in uniq_vals:
    #         val,ti,o = find_best(comp_adjusted,nc,c,uniq_curr,locs,start=iters)
    #         #print(val,ti)
    #         #val = find_best(comp,nc,c,uniq_curr,cc,start=iters)
    #         #val = np.sum(Mnn[np.where(nc == c)][:,np.where(nc==c)])
    #         count = 0
    #         fc = -1
    #         l = []
    #         #print(np.where(nc==c)[0])
    #         for i in np.where(nc==c)[0]:
    #             #print(i)
    #             clust = uniq_clust[i]
    #             #print(clust)
    #             fc = clust
    #             x=best_pmi[np.where(temp==clust)[0][0]]
    #             count+=x
    #             l.append((x,i))
    #         #print(c,val,count)
            
    #         if len(uniq_vals) < 8:
    #             print(val,ti,o)
    #             print(count,l)
            
    #         if val > count:
    #             for i in np.where(nc==c)[0]:
    #                 #print(i)
    #                 clust = uniq_clust[i]
    #                 curr_clusters[np.where(temp == clust)] = fc
    #                 #print(curr_clusters)
    #                 best_pmi[np.where(temp == clust)] = val
    #     # if len(uniq_vals) == 3:
    #     #     exit()
    #     _,curr_clusters = np.unique(curr_clusters,return_inverse=True)
    #     uniq_clust = np.unique(curr_clusters)
    # print(curr_clusters)
    # print(gt)
    # print(calc_pmi(comp,curr_clusters))
    # print(calc_pmi(comp,gt))
    # print(normalized_mutual_info_score(curr_clusters,gt))
   
   
    #print(normalized_mutual_info_score(gt,curr_clusters))


print('PMI')
# name = 'polblogs'
# file = 'Graphs/polblogs/network.pkl'
# name = 'cora'
# file = 'Graphs/cora/network.pkl'
# name = 'LFR'
# file = 'Graphs/LFR/network.pkl'
name = 'airport_ww'
file = 'Graphs/airport_ww/network.pkl'
# comp = list(np.load('computed_airport.npy'))
# name = 'wiki-fields'
# file = 'Graphs/wiki-fields/network.pkl'
# name = 'entsoe'
# file = 'Graphs/entsoe/network.pkl'
# name = 'custom'
# file = 'Graphs/custom/network.pkl'
# data = pd.read_csv('Graphs/airport_ww/node_info.csv')
# cnt = np.where(data['Continent'].to_numpy() == 'North America')[0]

#run(file=file,name=name,subgraph_locs=cnt)
run(file=file,name=name)
##############################################
# flag is first argument for runner function so for pmi flag is 'p'
# times='macro' allows you to change averaging scheme for lmepmi and mac
# SAVES: predictions in folder Predictions/{name}/pmi,Predictions/{name}/ac ... etc
# SAVES: Intermediate matrix exponentials in Predictions/{name}
#SAVES: Matrix exponential as comp_{name} in working directory.
# comp = run('p',file=file,name=name,times='micro',precomp=None)
# np.save('computed_{}'.format(name),comp)

# #comp = 'd'
# print('AC')
# run('ac',file=file,name=name,times='micro',precomp=comp)
# print()
# print()
# print('LMEPMI')
# run('lp',file=file,name=name,times='micro',precomp=comp)
# print()
# print('MAC')
# run('ma',file=file,name=name,times='micro',precomp = comp)