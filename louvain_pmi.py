import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import factorization as f
import scipy as sp
import time
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
Gcn = Gcn.add_nodes_from(G.nodes)
Gnc = Gnc.add_nodes_from(G.nodes)
nz = P.nonzero()
print(P.shape)
print(len(nz[0]),len(nz[1]))
for i in range(len(nz[0])):
    v = nz[0][i]+1
    u = nz[1][i]+1
    Gn.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1],pc=P[v-1,u-1]/Pi[v-1])
    Gn.add_edge(u,u,pn=P[u-1,u-1]/(Pi[u-1]),pc=P[u-1,u-1]/Pi[u-1])
    Gn.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1],pc = P[v-1,v-1]/Pi[u-1])
    Gcc.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1],pc=P[v-1,u-1]/Pi[v-1])
    Gcc.add_edge(u,u,pn=P[u-1,u-1]/Pi[u-1],pc=P[u-1,u-1]/Pi[u-1])
    Gcc.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1],pc=P[v-1,v-1]/Pi[v-1])
    Gcn.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1])
    Gcn.add_edge(u,u,pn=P[u-1,u-1]/Pi[u-1])
    Gcn.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1])
    Gnc.add_edge(v,u,pn=P[v-1,u-1]/Pi[v-1])
    Gnc.add_edge(u,u,pn=P[u-1,u-1]/Pi[u-1])
    Gnc.add_edge(v,v,pn=P[v-1,v-1]/Pi[v-1])
for v in Gcc.nodes():
    Gcc[v]['contains'] ={v}
    Gcc[v]['pi'] = Pi[v-1]



# for v in sorted(Gc.nodes):
#     print(Gc[v])
# for i in range(P.shape[0]):
Graphs = []
c = np.arange(P.shape[0])+1
Pia = Pi
while next_PMI-PMI > 0.001 or next_PMI == -1:
    for v in sorted(Gc.nodes):
        best=v
        best_inc = 0
        for u in Gn[v]:
            if u == v:
                continue
            c = c[u-1]
            cv = c[v-1]
            inc = -np.log(Pia[v-1]*Gcc[cv][cv]['pi']-Pia[v-1]*Gcc[cv]['pi'])+
            np.log(1+((Pia[v-1]*(Gn[v][u]['pc']+Gn[v][v]['pn'])+ Gcc[c]['pi'])/(Gcc[c]['pi']*Gcc[c][c]['pc']))) -
            np.log(1+((2*Pi[v-1]*Gcc[c]['pi']+Pi[v-1]**2)/(Gcc[c]['pi']**2)))
            if inc>best_inc:
                best = u
                best_inc = inc
        if best_inc == 0:
            continue
        c = c[best-1]
        cv = c[v-1]
        Pcc_new = (Gcc[c]['pi']*(Gcc[c][c]['pc']+Gcn[c][v]['pn']) + Pi[v-1]*(Gn[v][u]['pc']+Gn[v][v]['pn']))/(Gcc[c]['pi']+Pi[v-1])
        pic_new = Gcc[c]['pi']+Pi[v-1]
        Pcvcv_new = (Gcc[cv]['pi']*Gcc[cv][cv]['pc']-Pia[v-1]*Gn[v][u]['pc']-Gcc[cv]['pi']*Gcn[cv][n]['pn'])/(Gcc[cv]['pi']-Pia[v-1])
        picv_new = Gcc[cv]['pi']-Pia[v-1]
        

        

    

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




