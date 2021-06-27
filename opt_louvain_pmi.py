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
    print("here")
    print(PMI)
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
            #print(Mcc[cv][cv])
            
            #print(cv,Gcc.nodes[cv]['pi'], Pi[v-1],decr,Mcc[cv][cv],Gcc.nodes[cv]['pi']**2)
            Pcvcv_new = (Mcc[cv][cv]-Mnc[v][cv]-Mcn[cv][v]+Mnn[v][v])
            #print(Pcvcv_new)
            decr += Pcvcv_new
                # nodes = []
                # for n,i in enumerate(clusters):
                #     if i == cv:
                #         print(n+1,Pi[n])
                #         nodes.append(n+1)
                # print(nodes)
            #print(cv)
            #s = time.time()
            inc_pmi = Mcn[:,v]+Mnc[v]+Mnn[v][v]+decr
            inc_pmi[cv] = 0
            # print(Mcn)
            # print(Mnc)
            # print(inc_pmi)
            best = np.argmax(inc_pmi)
            best_inc = inc_pmi[best]
            
            #s = time.time()
            # for u in range(Mnn.shape[1]):
            #     c = clusters[u]
            #     if c == clusters[v]:
            #         continue
            #     new_Mcc = ((Mcc[c][c]+Mcn[c][v]) + (Mnc[v][c]+Mnn[v][v]))
            #     inc  = decr-Mcc[c][c]+new_Mcc
            #     #print(inc)
            #     #print(v,u,inc,new_Mcc,Mcc[c][c],next_PMI+inc,next_PMI)
            #     if np.isnan(inc) or np.isinf(inc):
            #         print(decr)
            #         #print(Pia[v-1],Gcc.nodes[c]['pi'],1+((Pia[v-1]*(Mnc[u][c]+Mnn[v][v])+ Gcc.nodes[c]['pi'])/(Gcc.nodes[c]['pi']*Mcc[c][c])))
            #         print(inc)
            #         print(new_Mcc)
            #         print("NAN")
            #         exit()
            #     #print(inc)
            #     if inc>t_inc:
            #         b_t = clusters[u]
            #         t_inc = inc
            # print(best,b_t,clusters[v])
            # print(best_inc,t_inc)
            # sec_time = time.time()-s
            # print("numpy:",first_time,"loop:",sec_time)
            # print("best",best_inc,"curr:",next_PMI)
            # print("test",t_inc)
            # print(best,t_best)
            # exit()
            if best_inc <= tol or best == clusters[v]:
                    continue
            else:
                increased = True
            c = best
            cv = clusters[v]
            # if v == 8:
            #     print(c)
            #     print(Mcc[c][c],Mnn[9][9])
            #     print(Mcn[c][v],Mnn[9][8])
            #     print(Mnc[v][c],Mnn[8][9])
            #     print(Gcc.nodes[c]['pi'],Pia[8])
            Mcc_new = ((Mcc[c][c]+Mcn[c][v]) +(Mnc[v][c]+Mnn[v][v]))
            #pic_new = Gcc.nodes[c]['pi']+Pi[v-1]
    
            Pcvcv_new = (Mcc[cv][cv]-Mnc[v][cv]-Mcn[cv][v]+Mnn[v][v])
            #picv_new = Gcc.nodes[cv]['pi']-Pi[v-1]
            Pvcv_new = Mnc[v][cv]-Mnn[v][v]
            Pvc_new = Mnc[v][c] + Mnn[v][v]
            Pcv_new = (Mcn[c][v]+Mnn[v][v])
            Pcvv_new= (Mcn[cv][v]-Mnn[v][v])
            # look at this section TODO
            Pccv_new = Mcc[c][cv]+Mnc[v][cv]-Mnn[v][v]-Mcn[c][v]
            Pcvc_new = Mcc[cv][c]-Mnc[v][c]+Mcn[cv][v]-Mnn[v][v]
            # if Pcvcv_new == 0:
            #     check = np.allclose(Mcc[cv],Mnc[v]) and np.allclose(Mcc[:,cv],Mcn[:,v])
            #     print(np.allclose(Mcc[cv],Mnc[v]))
            #     print(np.allclose(Mcc[:,cv],Mcn[:,v]))
            #     if not check:
            #         exit()
                # Mcc[cv][cv]-=Mnc[v]
                # print(Mcc[cv][cv])
            Mcc[c]+=Mnc[v]
            Mcc[cv]-=Mnc[v]
            Mcc[:,c]+=Mcn[:,v]
            Mcc[:,cv]-=Mcn[:,v]

            Mcn[c]+=Mnn[v]
            Mcn[cv]-=Mnn[v]
            Mnc[:,c]+=Mnn[:,v]
            Mnc[:,cv]-=Mnn[:,v]
            c1 = np.zeros(Mnn.shape[0])
            np.add.at(c1,clusters,Mnn[v]) # what you would think c1[clusters] += Mnn[v] would do
            c2 = np.zeros(Mnn.shape[0])
            np.add.at(c2,clusters,Mnn[:,v])
            # c3 = np.zeros(Mnn.shape[0])
            # for i in range(Mnn.shape[0]):
            #     c3[clusters[i]] += Mnn[i][v]
            # print(np.allclose(c2,c3))
           
            Mcc[c][c] = Mcc_new
            Mcc[cv][cv]= Pcvcv_new
            Mcc[c][cv]=Pccv_new
            Mcc[cv][c]=Pcvc_new
            # updates = set()
            # for x in Gn[v]:
            #     if x == v:
            #         continue
            #     cx = clusters[x-1]
            #     Gcc.add_edge(c,cx)
            #     Gcc.add_edge(cx,c)
            #     if cx not in Mcc[c]:
            #         Mcc[c][cx] = 0
            #         Mcc[cx][c] = 0
                    
            #     if cx not in updates:
            #         Mcc[c][cx] = Mcc[c][cx]
            #         if picv_new > tol:
            #             # print()
            #             # print(Mnn[v],x,cx)
            #             # print(Mcc[cv])
            #             # print()
            #             Mcc[cv][cx] = Mcc[cv][cx]
            #         else:
            #             Mcc[cv][cx] = 0
            #             Mcn[cv][x] = 0
            #         updates.add(cx)
            #     Mcc[c][cx] += (Mnn[v][x])
            #     if picv_new > tol:
            #         Mcc[cv][cx]-=Mnn[v][x]
            #     Mcc[cx][c]+=(Mnn[x][v])
            #     Mcc[cx][cv]-=Mnn[x][v]
            #     if x not in Mcn[c]:
            #         Mcn[c][x] = 0
            #     Mcn[c][x] = (Mcn[c][x]+Mnn[v][x])
            #     if picv_new > tol:
            #         Mcn[cv][x] = (Mcn[cv][x]-Mnn[v][x])
            #     else:
            #         Mcn[cv][x] = 0
            #     if c not in Mnc[x]:
            #         Mnc[x][c] = 0
            #     Mnc[x][c]+=Mnn[x][v]
            #     Mnc[x][cv] -= Mnn[x][v]
                
            # Mcc[c][c] = Mcc_new
            # Mcc[cv][cv] = Pcvcv_new
            # Mcn[c][v] = Pcv_new
            # Mcn[cv][v] = Pcvv_new
            # print("New cluster",c,":")
            # print("c to c:",Mcc[c][c])
            # if c not in Mnc[v]:
            #     Mnc[v][c] = 0
            # Mnc[v][c] += Mnn[v][v]
            # Mnc[v][cv]-=Mnn[v][v]
            # Gcc.nodes[c]['pi'] = pic_new
            # Gcc.nodes[cv]['pi'] = picv_new
            clusters[v] = c
            #print("moved node",v,"from cluster",cv,"to cluster", c)
            next_PMI += best_inc
            # clust = [[] for i in range(Mnn.shape[0])]
            # test_PMI = 0
            # for c,i in enumerate(clusters):
            #     clust[i].append(c)
            #     test_PMI += Mcc[c][c] 
            # sec_test = 0
            # for c in clust:
            #     for i in c:
            #         for j in c:
            #             sec_test+=Mnn[i][j]
         
            # if next_PMI- sec_test > tol:
            #     print(next_PMI,test_PMI,sec_test)
            #     # np.set_printoptions(threshold=sys.maxsize)
            #     Mnc_test = np.zeros(Mnc.shape)
            #     Mcn_test = np.zeros(Mcn.shape)
            #     for c,cl in enumerate(clust):
            #         for j in range(Mnn.shape[0]):
            #             s = 0
            #             for i in cl:
            #                 s+=Mnn[i][j]
            #             Mcn_test[c][j] = s
            #     for i in range(Mnn.shape[0]):
            #         for c,cl in enumerate(clust):
            #             s = 0
            #             for j in cl:
            #                 s+=Mnn[i][j]
            #             Mnc_test[i][c] = s
                
            #     print(np.allclose(Mcn_test,Mcn))
            #     print(np.allclose(Mnc,Mnc_test))
            #print("Increase:",best_inc,next_PMI)
        e = time.time()
        count+=1
        print()
        print("Run:",count,"Time:",e-s,"PMI:",next_PMI)
        # np.set_printoptions(threshold=sys.maxsize)
        # print(clusters)
        print()
    return (increased,Mcc,Mcn,Mnc,clusters,next_PMI)



def find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi):
    Graphs = []
    PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
    clusters = np.arange(P.shape[0])
    hierch_clusters = []
    increased = True 
    next_clusters = np.zeros(P.shape[0])
    valid_locs = np.arange(P.shape[0])
    s = time.time()
    increased,n_Mnn,n_Mcn,n_Mnc,n_clusters,PMI = onestep(PMI,clusters,0.0001,Mcc,Mcn,Mnn,Mnc,Pi)
    
    e =  time.time()
    print(e-s)
    print(n_clusters)
    uniq_vals,next_clusters = np.unique(n_clusters,return_inverse=True)
    next_Mnn = n_Mnn[uniq_vals][:,uniq_vals]
    
    #np.set_printoptions(threshold=sys.maxsize)
    #print(np.diag(n_Mnn))
    # t_Mcc = np.zeros(n_Mnn.shape)
    # clust = [[] for i in range(Mnn.shape[0])]
    # for c,i in enumerate(clusters):
    #     clust[i].append(c)
    # for c,cl in enumerate(clust):
    #     for c2,cl2 in enumerate(clust):
    #         s = 0
    #         for i in cl:
    #             for j in cl2:
    #                 s+=Mnn[i][j]
    #         t_Mcc[c][c2] = s
    # Mnc_test = np.zeros(n_Mnc.shape)
    # Mcn_test = np.zeros(n_Mcn.shape)
    # for c,cl in enumerate(clust):
    #     for j in range(n_Mnn.shape[0]):
    #         s = 0
    #         for i in cl:
    #             s+=Mnn[i][j]
    #         Mcn_test[c][j] = s
    # for i in range(n_Mnn.shape[0]):
    #     for c,cl in enumerate(clust):
    #         s = 0
    #         for j in cl:
    #             s+=Mnn[i][j]
    #         Mnc_test[i][c] = s
    # print((t_Mcc-n_Mnn).nonzero())
    # print(np.allclose(np.diag(t_Mcc),np.diag(n_Mnn)))
    # print(np.allclose(t_Mcc,n_Mnn))
    # print(np.allclose(Mcn_test,n_Mcn))
    # print(np.allclose(Mnc_test,n_Mnc))


    #np.save('Predictions/predicted_communities',n_clusters)
    print(next_clusters)
    print(PMI)
    #exit()
    hierch_clusters.append(next_clusters)
    count = 1
    print(count)
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
    return final_cluster

#G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
#G_1 = nx.karate_club_graph()
G = nx.read_gpickle('Graphs/airport_ww/network.pkl')
# G = nx.Graph()
# G.add_node(1)
# G.add_node(2)
# G.add_node(3)
# G.add_node(4)
# G.add_node(5)
# G.add_node(6)
# G.add_edge(1,2)
# G.add_edge(3,1)
# G.add_edge(2,3)
# G.add_edge(4,5)
# G.add_edge(4,6)
# G.add_edge(5,6)
# G.add_edge(3,4)




#G.add_edge(2,3)
# G.add_edge(1,2)
# G.add_nodes_from([i+1 for i in G_1.nodes])
# G.add_edges_from([(a[0]+1,a[1]+1) for a in G_1.edges])
print(G.number_of_nodes())
print(f.degree_vector(G).shape)
print(G.nodes)


d = f.degree_vector(G)
vert = np.arange(len(d))
print(d)
D= np.diag(1 / f.degree_vector(G))
A = f.adjacency_matrix(G)
P = np.asarray(np.matmul(D, A))
# D = sp.sparse.diags(1 / d)
# #print(D)
# A = nx.adjacency_matrix(G,sorted(G.nodes))

# print(sp.sparse.isspmatrix_csr(A))
# print(sp.sparse.isspmatrix_dia(D))
#s = time.time()
# P = D.dot(A)
# end = time.time()
# print(end-s)
#print(f.standard_random_walk_transition_matrix(G))
Pi = d/(np.sum(d))
print(Pi)
a = np.arange(10)+1
y0 =[]
y1= []
y2 = []
times = 10**np.linspace(-5,5,100)

times = [1]
#P = P.todense()
P_orig = np.copy(P)
print(P)
predictions = []
time_taken = []
for t in times:
    s = time.time()
    P = sp.linalg.expm((P_orig-np.eye(P_orig.shape[0]))*t)
    # print(P)
    # print(P.diagonal())
    # print(Pi)

    PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))

    print("INITIIAL:",PMI)
    # if PMI == np.log(0):
    #     PMI = -np.sum(np.log(Pi**2))

    next_PMI = -1
    # Gn = nx.Graph()
    # Gcc = nx.Graph()
    # Gcn = nx.Graph()
    # Gnc = nx.Graph()
    # Gn.add_nodes_from(G.nodes)
    # Gcc.add_nodes_from(G.nodes)
    

    # Gcn = Gcn.add_nodes_from(G.nodes)
    # Gnc = Gnc.add_nodes_from(G.nodes)
    #nz = P.nonzero()
    Mnn = np.log(P)-np.log(Pi)
    Mcc = np.copy(Mnn)
    Mcn = np.copy(Mnn)
    Mnc = np.copy(Mnn)
    print("Run:",t)
    final_cluster = find_clusters(Mcc,Mnn,Mcn,Mnc,P,Pi)
    e = time.time()
    time_taken.append(e-s)
    map_clust = {}
    mapped_clusters = []
    count = 1
    for i in final_cluster:
        if i not in map_clust:
            map_clust[i] = count
            count+=1
        mapped_clusters.append(map_clust[i])
    predictions.append(mapped_clusters)
print(mapped_clusters)
np.save('Predictions/predicted_communities_{}'.format(len(times)),np.array(predictions).T)
print(times)
# color_map = []
# # for i in final_cluster:
# #     if i == 31:
# #         color_map.append('green')
# #     if i == 5:
# #         color_map.append('red')
# #     if i == 13:
# #         color_map.append('orange')
# for i in final_cluster:
#     if i == 17:
#         color_map.append('red')
#     elif i == 27:
#         color_map.append('orange')
#     else:
#         color_map.append('green')

# nx.draw_spring(G,node_color=color_map,with_labels=True)
# plt.show()
