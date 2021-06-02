import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time
import copy
# TODO: Change Pcc and Pcn to only use numerator. Then, divide when neccesary.
min_inc = 0.0001
def onestep(PMI,clusters,min_inc,Pcc,Pcn,Pnn,Pnc,Gn,Gcc,Pi):
    tol =1e-15
    Pia = Pi
    increased = False
    next_PMI = -1
    print("here")
    print(PMI)
    while next_PMI-PMI > min_inc or next_PMI == -1:
        if next_PMI == -1:
            next_PMI = PMI
        else:
            PMI = next_PMI
        for v in Gn.nodes:
            best=v
            best_inc = 0
            cv = clusters[v-1]
            decr =   -np.log(Gcc.nodes[cv]['pi']*Pcc[cv][cv])+np.log(Gcc.nodes[cv]['pi']**2)
            if Gcc.nodes[cv]['pi'] - Pia[v-1] > tol:
                Pcvcv_new = (Gcc.nodes[cv]['pi']*Pcc[cv][cv]-Pia[v-1]*Pnc[v][cv]-Gcc.nodes[cv]['pi']*Pcn[cv][v]+Pia[v-1]*Pnn[v][v])/(Gcc.nodes[cv]['pi']-Pia[v-1])
                decr+= np.log((Gcc.nodes[cv]['pi']-Pia[v-1])*(Pcvcv_new)) - np.log((Gcc.nodes[cv]['pi']-Pia[v-1])**2)
                nodes = []
                for n,i in enumerate(clusters):
                    if i == cv:
                        print(n+1,Pia[n])
                        nodes.append(n+1)
                print(Pcvcv_new,Pia[v-1],Gcc.nodes[cv]['pi'],Pnc[v][cv],Pcn[cv][v],Pnn[v][v])
                print(cv,v)
                print("PCC",Pcc[cv][cv])
                print("Actual Values:")
                npcc = 0
                npcv = 0
                npvc = 0
                n_pi = 0
                print(nodes)
                for i in nodes:
                    n_pi += Pia[i-1]
                    if v in Pnn[i]:
                        npcv+=Pia[i-1]*Pnn[i][v]
                        npvc += Pia[v-1]*Pnn[v][i]
                    for j in nodes:
                        if j in Pnn[i]:
                            npcc+=Pia[i-1]*Pnn[i][j]

                npcc/=n_pi
                npcv/=n_pi
                npvc/=Pia[v-1]
                print(Pcn[cv][v],npcv)
                print(npcc,n_pi,npcv,npvc)
                ex = False
                if abs(npcc - Pcc[cv][cv]) > tol:
                    print("PCC",npcc,Pcc[cv][cv])
                    ex = True
                if abs(n_pi - Gcc.nodes[cv]['pi'] ) > tol:
                    print("pi",n_pi,Gcc.nodes[cv]['pi'])
                    ex = True
                if abs(npcv-Pcn[cv][v]) > tol:
                    print("PCN",npcv,Pcn[cv][v])
                    ex = True
                if abs(npvc-Pnc[v][cv]) > tol:
                    print("PCN",npvc,Pnc[v][cv])
                    ex = True
                if ex:
                    exit()
                # print((Pia[v-1]*(Pnn[v][v]+Pnn[v][8]))/(Pia[v-1]))
                # #print((Pia[7]*Pnn[8][v]+Pia[v-1]*Pnn[v][v])/(Pia[v-1]+Pia[7]))
                # print("PCC",Pcc[cv][cv])
                # print(v)
                # npc = (Pia[7]*(Pnn[8][8]+Pnn[8][v])+Pia[v-1]*(Pnn[v][v]+Pnn[v][8]))/(Pia[v-1]+Pia[7])
                # npcc = (Gcc.nodes[cv]['pi']*npc-Pia[v-1]*Pnc[v][cv]-Gcc.nodes[cv]['pi']*Pcn[cv][v]+Pia[v-1]*Pnn[v][v])/(Gcc.nodes[cv]['pi']-Pia[v-1])
                # print("npcc",npcc)
                # print("actual",(Pnn[8][8]))
                # print("npc",npc)
                # print((Pia[7]*(Pnn[8][8]+Pnn[8][v])+Pia[v]*(Pnn[v][v]+Pnn[v][8]))/(Pia[v-1]+Pia[7]))
                # print(Pia[v-1]*Pnc[v][cv],Gcc.nodes[cv]['pi'])
            for u in Gn[v]:
                if u == v:
                    continue
                c = clusters[u-1]
                if Gcc.nodes[c]['pi'] < tol:
                    continue
                new_Pcc = (Gcc.nodes[c]['pi']*(Pcc[c][c]+Pcn[c][v]) + Pia[v-1]*(Pnc[v][c]+Pnn[v][v]))/(Gcc.nodes[c]['pi']+Pia[v-1])

                inc  = decr+\
                    np.log((Gcc.nodes[c]['pi']+Pia[v-1])*new_Pcc)-np.log((Gcc.nodes[c]['pi']+Pia[v-1])**2)
                
                # np.log(1+((Pia[v-1]*(Pnc[v][c]+Pnn[v][v])+ Gcc.nodes[c]['pi'])/(Gcc.nodes[c]['pi']*Pcc[c][c]))) -\
                # np.log(1+((2*Pia[v-1]*Gcc.nodes[c]['pi']+Pi[v-1]**2)/(Gcc.nodes[c]['pi']**2)))
                if np.isnan(inc):
                    print(decr)
                    print(Pia[v-1],Gcc.nodes[c]['pi'],1+((Pia[v-1]*(Pnc[u][c]+Pnn[v][v])+ Gcc.nodes[c]['pi'])/(Gcc.nodes[c]['pi']*Pcc[c][c])))
                    print(new_Pcc)
                    print("NAN")
                    exit()
                #print(inc)
                if inc>best_inc:
                    best = u
                    best_inc = inc
            if best_inc <= 0:
                    continue
            else:
                increased = True
    
            c = clusters[best-1]
            cv = clusters[v-1]
            # if v == 8:
            #     print(c)
            #     print(Pcc[c][c],Pnn[9][9])
            #     print(Pcn[c][v],Pnn[9][8])
            #     print(Pnc[v][c],Pnn[8][9])
            #     print(Gcc.nodes[c]['pi'],Pia[8])
            Pcc_new = (Gcc.nodes[c]['pi']*(Pcc[c][c]+Pcn[c][v]) + Pia[v-1]*(Pnc[v][c]+Pnn[v][v]))/(Gcc.nodes[c]['pi']+Pia[v-1])
            pic_new = Gcc.nodes[c]['pi']+Pia[v-1]
            Pcvcv_new = 0
            if Gcc.nodes[cv]['pi']-Pia[v-1] <= tol:
                Pcvcv_new = 0
            else:
                Pcvcv_new = (Gcc.nodes[cv]['pi']*Pcc[cv][cv]-Pia[v-1]*Pnc[v][cv]-Gcc.nodes[cv]['pi']*Pcn[cv][v]+Pia[v-1]*Pnn[v][v])/(Gcc.nodes[cv]['pi']-Pia[v-1])
            picv_new = Gcc.nodes[cv]['pi']-Pia[v-1]
            Pvcv_new = Pnc[v][cv]-Pnn[v][v]
            Pvc_new = Pnc[v][c] + Pnn[v][v]
            Pcv_new = (Gcc.nodes[c]['pi']*Pcn[c][v]+Pia[v-1]*Pnn[v][v])/(pic_new)
            Pcvv_new = 0
            if abs(picv_new-0)<= tol:
                Pcvv_new = 0
            else:
                Pcvv_new= (Gcc.nodes[cv]['pi']*Pcn[cv][v]-Pia[v-1]*Pnn[v][v])/(picv_new)
            nodes = []
            for n,i in enumerate(clusters):
                if i == c:
                    print(n+1,Pia[n])
                    nodes.append(n+1)
            nodes.append(v)
            npcc = 0
            npcv = 0
            npvc = 0
            n_pi = 0
            old_npcv = 0
            print(nodes)
            print(v,c)
            for i in nodes:
                n_pi += Pia[i-1]
                if v in Pnn[i]:
                    npcv+=Pia[i-1]*Pnn[i][v]
                    npvc += Pia[v-1]*Pnn[v][i]
                    if i != v:
                        old_npcv+=Pia[i-1]*Pnn[i][v]
                for j in nodes:
                    if j in Pnn[i]:
                        npcc+=Pia[i-1]*Pnn[i][j]
            npcc/=n_pi
            npvc/=Pia[v-1]
            npcv/= n_pi
            old_npcv/=(n_pi-Pia[v-1])
            if 22 in Pcn[5]:
                print("PCN 5,22:",Pcn[5][22])
                print(Pnn[17])
                print(Pnn[5],Pia[4])
            print("PCN",npcv,Pcv_new)
            print("OLD",old_npcv,Pcn[c][v])
            ex = False
            if abs(npcc - Pcc_new) > tol:
                print("PCC",npcc,Pcc_new)
                ex = True
            if abs(n_pi - pic_new ) > tol:
                print("pi",n_pi,pic_new)
                ex = True
            if abs(npcv-Pcv_new) > tol:
                print("PCN",npcv,Pcv_new)
                print("OLD",old_npcv,Pcn[c][v])
                ex = True
            if abs(npvc-Pvc_new) > tol:
                #print(npvc-Pcv_new)
                print("PNC",npvc,Pvc_new)
                ex = True
            if ex:
                exit()
            
            updates = set()
            for x in Gn[v]:
                if x == v:
                    continue
                cx = clusters[x-1]
                Gcc.add_edge(c,cx)
                Gcc.add_edge(cx,c)
                if cx not in Pcc[c]:
                    if c in Pcc[cx]:
                        print("WEIRD")
                    Pcc[c][cx] = 0
                    Pcc[cx][c] = 0
                    
                if cx not in updates:
                    Pcc[c][cx] = Gcc.nodes[c]['pi']*Pcc[c][cx]/(pic_new)
                    if picv_new > tol:
                        Pcc[cv][cx] = Gcc.nodes[cv]['pi']*Pcc[cv][cx]/(picv_new)
                    else:
                        Pcc[cv][cx] = 0
                        Pcn[cv][x] = 0
                    updates.add(cx)
                Pcc[c][cx] += (Pi[v-1]*Pnn[v][x])/(pic_new)
                if picv_new > tol:
                    Pcc[cv][cx]-=Pi[v-1]*Pnn[v][x]/(picv_new)
                Pcc[cx][c]+=(Pi[x-1]*Pnn[x][v])/(Gcc.nodes[cx]['pi'])
                Pcc[cx][cv]-=(Pi[x-1]*Pnn[x][v])/(Gcc.nodes[cx]['pi'])
                if x not in Pcn[c]:
                    Pcn[c][x] = 0
                Pcn[c][x] = (Gcc.nodes[c]['pi']*Pcn[c][x]+Pia[v-1]*Pnn[v][x])/(pic_new)
                if picv_new > tol:
                    Pcn[cv][x] = (Gcc.nodes[cv]['pi']*Pcn[cv][x]-Pia[v-1]*Pnn[v][x])/(picv_new)
                else:
                    Pcn[cv][x] = 0
                if c not in Pnc[x]:
                    Pnc[x][c] = 0
                Pnc[x][c]+=Pnn[x][v]
                Pnc[x][cv] -= Pnn[x][v]
                
            nodes = []

         
            Pcc[c][c] = Pcc_new
            Pcc[cv][cv] = Pcvcv_new
            Pcn[c][v] = Pcv_new
            Pcn[cv][v] = Pcvv_new
            print("New cluster",c,":")
            print("c to c:",Pcc[c][c])
            if c not in Pnc[v]:
                Pnc[v][c] = 0
            Pnc[v][c] += Pnn[v][v]
            Pnc[v][cv]-=Pnn[v][v]
            if v not in Pcn[c]:
                Pcn[c][v] = 0
            Pcn[c][v] = (pic_new*Pcn[c][v]+Pia[v-1]*Pnn[v][v])/(pic_new)
            if picv_new <= tol:
                Pcn[cv][v] = 0
            else:
                Pcn[cv][v] = (Gcc.nodes[cv]['pi']*Pcn[cv][x]-Pia[v-1]*Pnn[v][v])/(picv_new)
            Gcc.nodes[c]['pi'] = pic_new
            Gcc.nodes[cv]['pi'] = picv_new
            clusters[v-1] = c
            print("moved node",v,"to", c)
            next_PMI+=best_inc
            print("Increase:",next_PMI)
    return (increased,Gcc,Pcc,Pcn,Pnc,clusters,next_PMI)




#G = nx.read_gpickle('Graphs/netsci/netsci_Gc.pkl')
G_1 = nx.karate_club_graph()
G = nx.Graph()
G.add_nodes_from([i+1 for i in G_1.nodes])
G.add_edges_from([(a[0]+1,a[1]+1) for a in G_1.edges])
print(G.number_of_nodes())
print(f.degree_vector(G).shape)

d = f.degree_vector(G)
vert = np.arange(len(d))
D = sp.sparse.diags(1 / d)
#print(D)
A = nx.adjacency_matrix(G,sorted(G.nodes))

print(sp.sparse.isspmatrix_csr(A))
print(sp.sparse.isspmatrix_dia(D))
s = time.time()
P = D.dot(A)
end = time.time()
print(end-s)
#print(f.standard_random_walk_transition_matrix(G))
Pi = 2*d/(np.sum(d))
print(Pi)
a = np.arange(10)+1
y0 =[]
y1= []
y2 = []
t = 1
P = P**2
#P_discrete = P
#P = sp.linalg.expm(P-np.eye(P.shape[0])*t)
print(P.diagonal())
PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
if PMI == np.log(0):
    PMI == -np.sum(np.log(Pi**2))

next_PMI = -1
Gn = nx.Graph()
Gcc = nx.Graph()
Gcn = nx.Graph()
Gnc = nx.Graph()
Gn.add_nodes_from(G.nodes)
Gcc.add_nodes_from(G.nodes)

Pnc = {}
Pcn = {}
Pcc = {}
Pnn = {}
final_clusters = {}
for i in G.nodes:
    Pnc[i] = {}
    Pcn[i] = {}
    Pcc[i] = {}
    Pnn[i] = {}
    final_clusters[i] = [i]
    print(type(Gcc))
    Gcc.nodes[i]['pi'] = Pi[i-1]
   

# Gcn = Gcn.add_nodes_from(G.nodes)
# Gnc = Gnc.add_nodes_from(G.nodes)
nz = P.nonzero()
print(P.shape)
print(len(nz[0]),len(nz[1]))
s = time.time()
for i in range(len(nz[0])):
    v = nz[0][i]+1
    u = nz[1][i]+1
    Gn.add_edge(v,u,p=P[v-1,u-1]/Pi[v-1],pc=P[v-1,u-1]/Pi[v-1])
    Gn.add_edge(u,u,p=P[u-1,u-1]/(Pi[u-1]),pc=P[u-1,u-1]/Pi[u-1])
    Gn.add_edge(v,v,p=P[v-1,v-1]/Pi[v-1],pc = P[v-1,v-1]/Pi[u-1])
    Gcc.add_edge(v,u,p=P[v-1,u-1]/Pi[v-1],pc=P[v-1,u-1]/Pi[v-1])
    Gcc.add_edge(u,u,p=P[u-1,u-1]/Pi[u-1],pc=P[u-1,u-1]/Pi[u-1])
    Gcc.add_edge(v,v,p=P[v-1,v-1]/Pi[v-1],pc=P[v-1,v-1]/Pi[v-1])
    # Pnc[v][u] = P[v-1,u-1]/Pi[v-1]
    # Pnc [v][v] = P[v-1,u-1]/Pi[v-1]
    # Pnc[u][u] = P[u-1,u-1]/Pi[u-1]
    # Pcn[v][u] = P[v-1,u-1]/Pi[v-1]
    # Pcn [v][v] = P[v-1,u-1]/Pi[v-1]
    # Pcn[u][u] = P[u-1,u-1]/Pi[u-1]
    # Pnn[v][u] = P[v-1,u-1]/Pi[v-1]
    # Pnn [v][v] = P[v-1,u-1]/Pi[v-1]
    # Pnn[u][u] = P[u-1,u-1]/Pi[u-1]
    # Pcc[v][u] = P[v-1,u-1]/Pi[v-1]
    # Pcc [v][v] = P[v-1,u-1]/Pi[v-1]
    # Pcc[u][u] = P[u-1,u-1]/Pi[u-1]
    Pnc[v][u] = P[v-1,u-1]
    Pnc [v][v] = P[v-1,v-1]
    Pnc[u][u] = P[u-1,u-1]
    Pnc[u][v] = P[u-1,v-1]
    Pcn[v][u] = P[v-1,u-1]
    Pcn [v][v] = P[v-1,v-1]
    Pcn[u][u] = P[u-1,u-1]
    Pcn[u][v] = P[u-1,v-1]
    Pnn[v][u] = P[v-1,u-1]
    Pnn [v][v] = P[v-1,v-1]
    Pnn[u][u] = P[u-1,u-1]
    Pnn[u][v] = P[u-1,v-1]
    Pcc[v][u] = P[v-1,u-1]
    Pcc [v][v] = P[v-1,v-1]
    Pcc[u][u] = P[u-1,u-1]
    Pcc[u][v] = P[u-1,v-1]
    


# for v in sorted(Gc.nodes):
#     print(Gc[v])
# for i in range(P.shape[0]):

Graphs = []
clusters = np.arange(P.shape[0])+1
hierch_clusters = []
Pia = Pi
increased = True 
next_clusters = np.zeros(P.shape[0])
valid_locs = np.arange(P.shape[0])
e = time.time()
print(e-s)
increased,n_Gnn,n_Pnn,n_Pcn,n_Pnc,n_clusters,PMI = onestep(PMI,clusters,0.0001,Pcc,Pcn,Pnn,Pnc,Gn,Gcc,Pi)
next_Pnn = copy.deepcopy(n_Pnn)
next_clusters = copy.deepcopy(n_clusters)
next_Gnn = n_Gnn
hierch_clusters.append((next_clusters,valid_locs))
count = 1
print(count)
while increased:
    increased = False
    valid_locs = []
    Pnn_new= {}
    Pcc_new={}
    Pcn_new = {}
    Pnc_new = {}
    Gcc_new = nx.Graph()
    valid = set()
    for a in next_clusters:
        valid.add(a)
    for v in list(next_Gnn.nodes):
        if v not in valid:
            next_Gnn.remove_node(v)
            continue
        Pnn_new[v] = next_Pnn[v]
        Pnc_new[v] = next_Pnn[v]
        Pcn_new[v] = next_Pnn[v]
        Pcc_new[v] = next_Pnn[v]
        valid_locs.append(v-1)
        Gcc_new.add_node(v)
        Gcc_new.add_edges_from([(v,u) for u in next_Gnn[v]])
        Gcc_new.nodes[v]['pi'] = next_Gnn.nodes[v]['pi']
        Pi[v-1] = next_Gnn.nodes[v]['pi']
        increased,next_Gnn,next_Pnn,next_Pcn,next_Pnc,next_clusters,PMI = onestep(PMI,clusters,0.0001,Pcc_new,Pcn_new,Pnn_new,Pnc_new,next_Gnn,Gcc_new,Pi)
        count+=1
        print(count)
final_cluster = np.zeros(P.shape[0])
print(hierch_clusters)
for v in Gn.nodes:
    c = v
    for p in hierch_clusters:
        c = p[0][c-1]
    final_cluster[v-1] = c
print(final_cluster)
    
    


    





    
    







        

        

    

# while  next_PMI - PMI > 0.001 or next_PMI == -1:
#     pass


# Pd = P.toarray()
# for i in a:
#     s = time.time()
#     P**i
#     e = time.time()
#     y0.append(e-s)
#     s = time.time()
#     # P.matrix_power(i)
#     e = time.time()
#     y1.append(e-s)
#     s = time.time()
#     np.linalg.matrix_power(Pd,i)
#     e = time.time()
#     y2.append(e-s)

# fig,axs = plt.subplots()
# axs.scatter(a,y0,color='red')
# axs.scatter(a,y2,color='blue')
# #axs[1].plot(a,y1)
# #axs[0].scatter(a,y2)
# plt.show()



    









# print(A)
# print()
# print()
# print(Ap)
# print(Ap.todense())




