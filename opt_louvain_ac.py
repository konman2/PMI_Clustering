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
            #     if u == v:
            #         continue
            #     c = clusters[u]
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
            #     if inc>best_inc:
            #         best = u
            #         best_inc = inc
            # sec_time = time.time()-s
            # print("numpy:",first_time,"loop:",sec_time)
            # print("best",best_inc,"curr:",next_PMI)
            # print("test",t_inc)
            # print(best,t_best)
            # exit()
            if best_inc <= tol or clusters[best] == clusters[v]:
                    continue
            else:
                increased = True
            c = clusters[best]
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
            Mcn[c]+=Mnn[v]
            Mcn[cv]-=Mnn[v]
            Mnc[:,c]+=Mnn[:,v]
            Mnc[:,cv]-=Mnn[:,v]
            Mcc[c][c] = Mcc_new
            Mcc[cv][cv]= Pcvcv_new
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
            next_PMI+=best_inc
            #print("Increase:",best_inc,next_PMI)
        e = time.time()
        count+=1
        print()
        print("Run:",count,"Time:",e-s,"PMI:",next_PMI)
        #np.set_printoptions(threshold=sys.maxsize)
        #print(clusters)
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
    np.save('Predictions/predicted_communities',n_clusters)
    exit()
    next_Mnn = copy.deepcopy(n_Mnn)
    next_clusters = copy.deepcopy(n_clusters)
    print(next_clusters)
    print(PMI)
    #exit()
    next_Gnn = n_Gnn
    hierch_clusters.append((next_clusters,valid_locs))
    count = 1
    print(count)
    while increased:
        increased = False
        valid_locs = []
        Mnn_new= {}
        Mcc_new={}
        Mcn_new = {}
        Mnc_new = {}
        Gcc_new = nx.Graph()
        valid = set()
        for a in next_clusters:
            valid.add(a)
        for v in list(next_Gnn.nodes):
            if v not in valid:
                next_Gnn.remove_node(v)
                continue
            Mnn_new[v] = next_Mnn[v]
            Mnc_new[v] = next_Mnn[v]
            Mcn_new[v] = next_Mnn[v]
            Mcc_new[v] = next_Mnn[v]
            
            valid_locs.append(v-1)
            Gcc_new.add_node(v)
            Gcc_new.add_edges_from([(v,u) for u in next_Gnn[v]])
            Gcc_new.nodes[v]['pi'] = next_Gnn.nodes[v]['pi']
            Pi[v-1] = next_Gnn.nodes[v]['pi']
        
        increased,next_Gnn,next_Mnn,next_Mcn,next_Mnc,next_clusters,PMI = onestep(PMI,clusters,0.0001,Mcc_new,Mcn_new,Mnn_new,Mnc_new,next_Gnn,Gcc_new,Pi)
        count+=1
        print("Number of runs:", count)
        
    final_cluster = np.zeros(P.shape[0])
    #print(hierch_clusters)
    for v in Gn.nodes:
        c = v
        for p in hierch_clusters:
            c = p[0][c-1]
        final_cluster[v-1] = c
    #print(final_cluster)
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
    Gn = nx.Graph()
    Gcc = nx.Graph()
    Gcn = nx.Graph()
    Gnc = nx.Graph()
    Gn.add_nodes_from(G.nodes)
    Gcc.add_nodes_from(G.nodes)
    

    # Gcn = Gcn.add_nodes_from(G.nodes)
    # Gnc = Gnc.add_nodes_from(G.nodes)
    nz = P.nonzero()
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
        mapped_clusters = map_clust[i]
    predictions.append(mapped_clusters)

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
