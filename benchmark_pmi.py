import numpy as np
import matplotlib.pyplot as plt

taus = 50
filename = f'Predictions/airport/pmi/predicted_communities_{taus}.npy'
filename_ac = f'Predictions/airport/ac/predicted_communities_{taus}.npy'

filename_lmepmi = f'Predictions/airport/lmepmi/predicted_communities_{taus}.npy'
filename_mac = f'Predictions/airport/mac/predicted_communities_{taus}.npy'

times = np.load('Predictions/airport/pmi/times.npy')
matr = np.load('computed_airport.npy')
macro_comms = np.load('Graphs/airport_ww/macro_comms.npy')
micro_comms = np.load('Graphs/airport_ww/micro_comms.npy')
pmi_vals = np.load('Predictions/airport/pmi/values.npy')
ac_vals = np.load('Predictions/airport/ac/values.npy')
Pi = np.load('Pis.npy')

predictions = np.load(filename)
predictions_ac = np.load(filename_ac)
def calculate_pmi(matr,times,Pi,pred):
    vals = []
    for it,t in enumerate(times):
        P = matr[it]
        #pred = predictions[:,it]
        Mnn = np.log(P)-np.log(Pi)
        #num_clusters = len(np.unique(micro_comms))
        pmi = 0
        for v in range(Mnn.shape[0]):
            c = pred[v]
            pmi+=np.sum(Mnn[v][np.nonzero(pred == c)])
        vals.append(pmi)
        #print(pmi,vals[it])
        #print(micro_pmi)
    return np.array(vals)
def calculate_ac(matr,times,Pi,pred):
    vals = []
    for it,t in enumerate(times):
        P = matr[it]
        #pred = predictions[:,it]
        Mnn = ((P-Pi).T*Pi).T
        #num_clusters = len(np.unique(micro_comms))
        ac = 0
        for v in range(Mnn.shape[0]):
            c = pred[v]
            ac+=np.sum(Mnn[v][np.nonzero(pred == c)])
        vals.append(ac)
        #print(ac,ac_vals[it])
        #print(pmi,vals[it])
        #print(micro_pmi)
    return np.array(vals)

#calculate_ac(matr,times,Pi,predictions_ac)
#exit()
pmis = calculate_pmi(matr,times,Pi,micro_comms)
# print(pmis)
# print(pmi_vals)
# exit()
acs = calculate_ac(matr,times,Pi,micro_comms)
print(acs)
print(ac_vals)
plt.plot(np.log10(times),pmis,label='micro community')
plt.plot(np.log10(times),pmi_vals,label='predicted community')
plt.legend()
plt.figure()
plt.plot(np.log10(times),acs,label='micro community')
plt.plot(np.log10(times),ac_vals,label='predicted community')
plt.legend()
plt.show()






