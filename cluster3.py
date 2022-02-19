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

min_inc = 0.0001
def onestep(PMI,clusters,min_inc,Mcc,Mcn,Mnn,Mnc,Pi):
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

def find_best2(locs,comp,start=0):
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
    #return (best/np.sum(locs),time,locs)
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

    # low = start
    # hi = len(comp)-1
    # best = -1
    # while(hi>low):
    #     mid = (hi+low)//2
    #     print(mid,hi,low)
    #     v1 = func(comp[mid+1],locs,nc,c,uniq_curr)
    #     v2 = func(comp[mid],locs,nc,c,uniq_curr)
    #     if v2 < v1:
    #         low = mid+1
    #     else:
    #         hi = mid
    #     print(v1,v2)
    #     best = max(best,max(v2,v1))
    # print(best,func(comp[0],locs,nc,c,uniq_curr))
    # return best

    # M = comp[start]
    # Mnn = np.zeros((len(uniq_curr),len(uniq_curr)))
    # for i in uniq_curr:
    #     for j in uniq_curr:
    #         Mnn[i][j] = np.sum(M[np.where(cc == i)][:,np.where(cc==j)])
    # best = np.sum(Mnn[np.where(nc == c)][:,np.where(nc==c)])
    # for t in range(start+1,len(comp)):
    #     M = comp[t]
    #     Mnn = np.zeros((len(uniq_curr),len(uniq_curr)))
    #     for i in uniq_curr:
    #         for j in uniq_curr:
    #             Mnn[i][j] = np.sum(M[np.where(cc == i)][:,np.where(cc==j)])
    #     val = np.sum(Mnn[np.where(nc == c)][:,np.where(nc==c)])
    #     if val>best:
    #         best = val
    # return best

def calc_pmi(comp,clust,locs=None):
    if str(type(locs)) == '<class \'NoneType\'>':
        locs = np.ones(len(clust), dtype=bool)
    uniq_curr,cc = np.unique(clust,return_inverse=True)
    total = 0
    clust_scores = []
    for c in uniq_curr:
        Mnn = comp[0]
        clust_locs = (clust==c) & locs
        best = np.sum(Mnn[clust_locs][:,clust_locs])
        best_time = 0
        for t in range(1,len(comp)):
            Mnn = comp[t]
            val = np.sum(Mnn[clust_locs][:,clust_locs])
            if val > best:
                best = val
                best_time = t
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

def run(flag='p',file='Graphs/airport_ww/network.pkl',name='airport',times='micro',precomp=None):
    #G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
    #G = nx.karate_club_graph()
    # gt_micro = []
    # gt_macro = []
    
    if name == 'airport':
        name = 'airport_ww'
    gt_micro = np.load(f'Graphs/{name}/micro_comms.npy')
    gt_macro = np.load(f'Graphs/{name}/macro_comms.npy')


    G = nx.read_gpickle(file)
    # gt = np.load(f'Graphs/{name}/micro_comms.npy')
    # assert len(gt) == len(G.nodes)

    # G,gt_micro,gt_macro = arlei_graph()
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

    sigma = 0
    d = f.degree_vector(G)
    vert = np.arange(len(d))
    #print(d)
    # D= np.diag(1 / f.degree_vector(G))
    # A = f.adjacency_matrix(G)
    #P = np.asarray(np.matmul(D, A))
    #P = f.pagerank_transition_matrix(G,mu=0.1)
    #print(np.nonzero(P==0))
    P = f.standard_random_walk_transition_matrix(G)
    Pi = d/(np.sum(d))
    #pi = f.stationary_distribution(P)
    #Pi = np.diag(pi)
    #print(Pi)
    #np.save('Pis_{}'.format(name),Pi) 
    a = np.arange(10)+1
    y0 =[]
    y1= []
    y2 = []
   
    #times = 10**np.linspace(-2,2,100)

    #times = np.linspace(10**(-3),10**3,2000)
    times = np.linspace(10**(-2),10**2,200)
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
        #print(t)
        P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
        M = np.log(P+sigma)-np.log(Pi+sigma)
        comp.append(M)
    comp = np.array(comp)
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


    for iters,t in enumerate(times):
        M = comp[iters]
        uniq_curr,cc = np.unique(curr_clusters,return_inverse=True)
        uniq_clust = uniq_curr
        Mnn = np.copy(M)
        #Mnn = np.zeros((len(uniq_curr),len(uniq_curr)))
        #locs = np.array([cc==i for i in uniq_curr])
        #comp_adjusted  = []
        #s1 = time.time()
        # for it in range(len(times)):
        #     for i in uniq_curr:
        #         # Mnn[i][cc[(cc==i)]]
        #         for j in uniq_curr:
        #             Mnn[i][j] = np.sum(M[locs[i]][:,locs[j]])
        #comp_adjusted = np.zeros((comp.shape[0],len(uniq_curr),len(uniq_curr)))
        # for i in uniq_curr:
        #         # Mnn[i][cc[(cc==i)]]
        #     for j in uniq_curr:
        #         Mnn[i][j] = np.sum(M[locs[i]][:,locs[j]])
            #comp_adjusted.append(Mnn) 
        #e1 = time.time()
        # for i in uniq_curr:
        #     for j in uniq_curr:
        #         # if iters > 0:
        #         #     print(iters)
        #         #     print(comp_adjusted[:,i,j].shape)
        #         #     print(np.sum(comp[:,locs[i]][:,:,locs[j]],axis=1).shape)
        #         #     print(np.where(locs[i]),np.where(locs[j]))
        #         #     print(comp[:,locs[i]][:,:,locs[j]])
        #         comp_adjusted[:,i,j]=np.sum(comp[:,locs[i]][:,:,locs[j]],axis=(1,2))
        # print(comp_adjusted[iters])
        # print()
        # print(Mnn)
        # print()
        #print(Mnn-comp_adjusted[iters])
        #continue
        #Mnn = comp_adjusted[iters]
        #Mnn = comp_adjusted[iters]
        #print((Mnn == comp[0]))
        
        #s2 = time.time()

        
        Mcc = np.copy(Mnn)
        Mcn = np.copy(Mnn)
        Mnc = np.copy(Mnn)
        final_cluster,pmi = find_clusters(Mcc,Mnn,Mcn,Mnc,Pi,np.sum(np.diag(Mnn)))
        uniq_vals,nc = np.unique(final_cluster,return_inverse=True)
        print(nc,len(uniq_vals), t, file = sourceFile)
        predictions.append(nc)
        
    off_limits = set()
    signal = np.zeros(P.shape[0])
    count = 0
    sourceFile.close()
    print("started")


    # loc = predictions[299] == 0
    
    # assign = predictions[299][loc]
    # assign+=P.shape[0]+count
    # for i in range(len(predictions)):
    #         predictions[i][loc] = assign
    # signal[loc] = 1
    # count = 1
    # while(np.sum(signal)<P.shape[0]):
    #     loc_s = set()
    #     loc_c = []
    #     loc_val = []
    #     loc_list = []
    #     #print(loc_s)
    #     for i, clustering in enumerate(reversed(predictions)):
    #         clusters = np.unique(clustering)
    #         #print(loc_s)
    #         # if count > 0:
    #         #     print(clusters)
    #         for c in clusters:
    #             if c<P.shape[0]:
    #                 l = (clustering == c)
    #                 # if count > 0:
    #                 #     print(c)
    #                 #     print(clustering)
    #                 #     print(l)
    #                 uniq_sum = str(l)
    #                 if uniq_sum not in loc_s:
    #                     loc_s.add(uniq_sum)
    #                     #print(l)
    #                     subset = False
    #                     to_remove = []
    #                     for k,loc in enumerate(loc_list):
    #                         if np.all(loc & l == l):
    #                             subset = True
    #                             break
    #                         if np.all(loc&l == loc):
    #                             to_remove.append(k)
    #                     if not subset:
    #                         for k in to_remove:
    #                             loc_list.remove(k)
    #                         loc_list.append(l)
                           
    #     for l in loc_list:
    #         if np.sum(l) == 1:
    #             loc_val.append(np.sum(comp[0][l][:,l]))
    #         else:
    #             loc_val.append(0)
    #         loc_c.append(0)

    #         # if count > 1:
    #         #      print(np.sum(loc_list,axis=1))
    #     loc_val = np.array(loc_val)
    #     #print(len(loc_list))
        
    #     for i,loc in enumerate(loc_list):
    #         if loc_val[i] != 0:
    #             print(loc_val[i],loc_c[i])
    #             continue
    #         best = 0
    #         best_c = 0
    #         for ind,clustering in enumerate(predictions):
    #             #print(np.where(loc)[0])
    #             val,_ = calc_pmi(comp,clustering,loc)
    #             val/=(np.sum(loc))
    #             #print()
    #             if val > best:
    #                 best = val
    #                 best_c = ind
    #         print(best,best_c,i)
    #         loc_val[i] = best
    #         loc_c[i] = best_c
    #     b_ind = np.argmax(loc_val)
    #     clust = np.copy(predictions[loc_c[b_ind]])
    #     loc = loc_list[b_ind]
    #     uniq_vals,assign = np.unique(clust[loc],return_inverse=True)
    #     assign+=(P.shape[0]+count)
    #     count+=len(uniq_vals)
    #     for i in range(len(predictions)):
    #         predictions[i][loc] = assign
    #     #print(predictions[395])
    #     signal[loc] = 1
    #     print(np.sum(signal))
    # #print(predictions)
    # _,o = np.unique(predictions[0],return_inverse=True)
    # print(o)
    # print(calc_pmi(comp,o))
    # print(gt)
    # print(calc_pmi(comp,gt))
    predictions = np.array(predictions)
    total = 0
    while(total<P.shape[0]):
        best_l = -1
        best_ind  = -1
        best = 0
        #proc = set()
        #valid_locs = np.array([])
        valid_locs = []
        # s = time.time()
        for i,clustering in enumerate(predictions):
            clusters = np.unique(clustering)
            #print(clustering)
            #print
            for c in clusters:
                if c >= P.shape[0]:
                    continue
                l = clustering == c
                #z = str(l)
                l = l.reshape(1,P.shape[0])
                #print(len(valid_locs))
                if len(valid_locs) == 0:
                    valid_locs = np.copy(l)
                elif np.max(np.sum(valid_locs == l,axis=1)) < P.shape[0]:
                    valid_locs = np.append(valid_locs,l,axis=0)
                

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
        vals = np.apply_along_axis(find_best2,1,valid_locs,comp)
        
        #print(vals)
        # exit()
        best_i = np.argmax(vals)
        best = vals[best_i]
        # assert np.sum(best_l == valid_locs[best_i]) == P.shape[0]
        best_l = valid_locs[best_i]
        # e = time.time()
        # print(e-s)
        s = time.time()
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
        print(best,np.sum(best_l))
        #print(np.sum(best_l))
        # e = time.time()
        # print(e-s)
        # print()
        #print(np.sum(signal))
        #print(best,best_l)
        #print((np.sum(l),i) for i,l in enumerate(locs))
        #print(predictions[0])
    _,o = np.unique(predictions[0],return_inverse=True)
    outfile = open('results.txt', 'w')
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
    outfile.close()
  



            


        
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