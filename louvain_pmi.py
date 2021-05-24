import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time

min_inc = 0.0001
def onestep(PMI,clusters,min_inc,Pcc,Pcn,Pnn,Pnc,Gn,Gcc,Pi):
    tol =1e-15
    Pia = Pi
    increased = False
    while next_PMI-PMI > min_inc or next_PMI == -1:
        PMI = next_PMI
        for v in sorted(Gn.nodes):
            best=v
            best_inc = 0
            cv = clusters[v-1]
            decr =   -np.log(Gcc[c]['pi']*Pcc[cv][cv])+np.log(Gcc[c]['pi']**2)
            if Gcc[cv]['pi'] - Pia[v-1] > tol:
                Pcvcv_new = (Gcc[cv]['pi']*Pcc[cv][cv]-Pia[v-1]*Pnc[v][cv]-Gcc[cv]['pi']*Pcn[cv][v])/(Gcc[cv]['pi']-Pia[v-1])
                decr+= np.log((Gcc[cv]['pi']-Pia[v-1])*(Pcvcv_new)) - np.log((Gcc[cv]['pi']-Pia[v-1])**2)
            
            for u in Gn[v]:
                if u == v:
                    continue
                c = clusters[u-1]
                inc  = decr+\
                np.log(1+((Pia[v-1]*(Pnc[u][c]+Pnn[v][v])+ Gcc[c]['pi'])/(Gcc[c]['pi']*Pcc[c][c]))) -\
                np.log(1+((2*Pi[v-1]*Gcc[c]['pi']+Pi[v-1]**2)/(Gcc[c]['pi']**2)))
                if inc>best_inc:
                    best = u
                    best_inc = inc
            if best_inc <= 0:
                    continue
            else:
                increased = True
    
            c = clusters[best-1]
            cv = clusters[v-1]
            Pcc_new = (Gcc[c]['pi']*(Pcc[c][c]+Pcn[c][v]) + Pi[v-1]*(Pnc[v][c]+Pnn[v][v]))/(Gcc[c]['pi']+Pia[v-1])
            pic_new = Gcc[c]['pi']+Pia[v-1]
            Pcvcv_new = 0
            if Gcc[cv]['pi']-Pia[v-1] <= tol:
                Pcvcv_new = 0
            else:
                Pcvcv_new = (Gcc[cv]['pi']*Pcc[cv][cv]-Pia[v-1]*Pnc[v][cv]-Gcc[cv]['pi']*Pcn[cv][v])/(Gcc[cv]['pi']-Pia[v-1])
            picv_new = Gcc[cv]['pi']-Pia[v-1]
            Pvcv_new = Pnc[v][cv]-Pnn[v][v]
            Pvc_new = Pnc[v][c] +Pnn[v][v]
            Pcv_new = (Gcc[c]['pi']*Pcn[c][v]+Pia[v-1]*Pnn[v][v])/(pic_new)
            if abs(picv_new-0)<= tol:
                Pcvv_new = 0
            else:
                Pcvv_new= (Gcc[cv]['pi']*Pcn[cv][v]-Pia[v-1]*Pnn[v][v])/(picv_new)
            updates = set()
            for x in Gn[v]:
                if x == v:
                    continue
                cx = clusters[x-1]
                Gcc.add_edge(c,cx)
                Gcc.add_edge(cx,c)
                
                if cx not in Pcc[c]:
                    Pcc[c][cx] = 0
                    Pcc[cx][c] = 0
                    
                if cx not in updates:
                    if cx not in Pcc:
                        Pcc[c][cx] = 0
                    Pcc[c][cx] = Gcc[c]['pi']*Pcc[c][cx]/(pic_new)
                    if picv_new <= tol:
                        Pcc[cv][cx] = Gcc[cv]['pi']*Pcc[cv][cx]/(picv_new)
                    else:
                        Pcc[cv][cx] = 0
                        Pcn[cv][x] = 0
                    updates.add(cx)
                Pcc[c][cx] += Pi[v-1]*Pnn[v][x]/(pic_new)
                if picv_new > tol:
                    Pcc[cv][cx]-=Pi[v-1]*Pnn[v][x]/(picv_new)
                Pcc[cx][c]+=Pi[x-1]*Pnn[x][v]/(Gcc[cx]['pi'])
                Pcc[cx][cv]-=Pi[x-1]*Pnn[x][v]/(Gcc[cx]['pi'])
                if x not in Pcn[c]:
                    P[c][x] = 0
                Pcn[c][x] = (Gcc[c]['pi']*Pcn[c][x]+Pia[v-1]*Pnn[v][x])/(pic_new)
                if picv_new > tol:
                    Pcn[cv][x] = (Gcc[cv]['pi']*Pcn[cv][x]-Pia[v-1]*Pnn[v][x])/(picv_new)
                if c not in Pnc[x]:
                    Pnc[x][c] = 0
                Pnc[x][c]+=Pnn[x][v]
                Pnc[x][cv] -= Pnn[x][v]
            Pcc[c][c] = Pcc_new
            Pcc[cv][cv] = Pcvcv_new
        
            if c not in Pnc[v][c]:
                Pnc[v][c] = 0
            Pnc[v][c] += Pnn[v][v]
            Pnc[v][cv]-=Pnn[v][v]
            if v not in Pcn[c]:
                Pcn[c][v] = 0
            Pcn[c][v] = (Gcc[c]['pi']*Pcn[c][v]+Pia[v-1]*Pnn[v][v])/(pic_new)
            if picv_new <= tol:
                Pcv[cv][v] = 0
            else:
                Pcn[cv][v] = (Gcc[cv]['pi']*Pcn[cv][x]-Pia[v-1]*Pnn[v][v])/(picv_new)
            Gcc[c]['pi'] = pic_new
            Gcc[cv]['pi'] = picv_new
            clusters[v-1] = c
            next_PMI+=best_inc
    return (increased,Gcc,Pcc,Pcn,Pnc,clusters,next_PMI)




G = nx.read_gpickle('netsci/netsci_Gc.pkl')
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
P = P**2
print(P.diagonal())
PMI = np.sum(np.log(Pi*P.diagonal()))-np.sum(np.log(Pi**2))
if PMI == np.log(0):
    PMI == -np.sum(np.log(Pi**2))

next_PMI = -1
Gn = nx.Graph()
Gcc = nx.Graph()
Gcn = nx.Graph()
Gnc = nx.Graph()
Gn = Gn.add_nodes_from(G.nodes)
Gcc = Gcc.add_nodes_from(G.nodes)
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
    final_clusters[i] = {i}

# Gcn = Gcn.add_nodes_from(G.nodes)
# Gnc = Gnc.add_nodes_from(G.nodes)
nz = P.nonzero()
print(P.shape)
print(len(nz[0]),len(nz[1]))
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
    Pnc [v][v] = P[v-1,u-1]
    Pnc[u][u] = P[u-1,u-1]
    Pcn[v][u] = P[v-1,u-1]
    Pcn [v][v] = P[v-1,u-1]
    Pcn[u][u] = P[u-1,u-1]
    Pnn[v][u] = P[v-1,u-1]
    Pnn [v][v] = P[v-1,u-1]
    Pnn[u][u] = P[u-1,u-1]
    Pcc[v][u] = P[v-1,u-1]
    Pcc [v][v] = P[v-1,u-1]
    Pcc[u][u] = P[u-1,u-1]
    

    




    # Gcn.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1])
    # Gcn.add_edge(u,u,pn=P[u-1,u-1]/Pi[u-1])
    # Gcn.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1])
    # Gnc.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1])
    # Gnc.add_edge(u,u,pn=P[u-1,u-1]/Pi[u-1])
    # Gnc.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1])
for v in Gcc.nodes():
    Gcc[v]['pi'] = Pi[v-1]





# for v in sorted(Gc.nodes):
#     print(Gc[v])
# for i in range(P.shape[0]):
Graphs = []
clusters = np.arange(P.shape[0])+1
hierch_clusters = []
Pia = Pi
increased = True 
next_clusters = np.zeros(P.shape[0])
while increased:
    
    increased = False
    increased,next_Gnn,next_Pnn,next_Pcn,next_Pnc,next_clusters,PMI = onestep(PMI,clusters,0.0001,Pcc,Pcn,Pnn,Pnc,Gn,Gcc,Pi)
    hierch_clusters.append(next_clusters)
    Gnn_new = nx.Graph()
    Pnn_new= {}
    Pcc={}
    Pcn = {}
    Pnn={}
    for v in next_Gnn.nodes:
        if clusters[v-1] not in final_clusters:
            next_clusters[clusters[v-1]] = []
        next_clusters[c[v-1]]+=final_clusters[v]
        c = clusters[v-1]
        next_clusters[c-1]=c
        Gnn_new.add_node(c)
        for u in Gnn[c]:
            pass



    

    final_clusters = next_clusters

    
    


    





    
    







        

        

    

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




