import numpy as np
import matplotlib.pyplot as plt
name = 'entsoe'
taus = 50
filename = f'Predictions/{name}/pmi/predicted_communities_{taus}.npy'
filename_ac = f'Predictions/{name}/ac/predicted_communities_{taus}.npy'

filename_lmepmi = f'Predictions/{name}/lmepmi/predicted_communities_{taus}.npy'
filename_mac = f'Predictions/{name}/mac/predicted_communities_{taus}.npy'

times = np.load(f'Predictions/{name}/pmi/times.npy')
matr = np.load(f'computed_{name}.npy')
pmi_vals = np.load(f'Predictions/{name}/pmi/values.npy')
ac_vals = np.load(f'Predictions/{name}/ac/values.npy')
Pi = np.load(f'Pis_{name}.npy')
name2 = name
if name == 'airport':
    name2 = 'airport_ww'
macro_comms = np.load(f'Graphs/{name2}/macro_comms.npy')
micro_comms = np.load(f'Graphs/{name2}/micro_comms.npy')


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
        print(it,pmi)
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
        print(ac,it)
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

tgt = np.log10(times[np.argmax(pmis)])
tp = np.log10(times[np.argmax(pmi_vals)])
tgt_ac = np.log10(times[np.argmax(acs)])
tp_ac = np.log10(times[np.argmax(ac_vals)])
plt.title('PMI vs Markov Time for Ground Truth and predicted communities')
plt.xlabel('Markov Time')
plt.ylabel('PMI')
plt.plot(np.log10(times),pmis,label=f'micro community, best at {tgt}')
plt.plot(np.log10(times),pmi_vals,label=f'predicted community, best at {tp}')
plt.legend()
#plt.savefig(f'fig/PMI_{name}_gt')
plt.figure()
plt.title('AC vs Markov Time for Ground Truth and predicted communities')
plt.xlabel('Markov Time')
plt.ylabel('AC')
plt.plot(np.log10(times),acs,label=f'micro community, best at {tgt_ac}')
plt.plot(np.log10(times),ac_vals,label=f'predicted community, best at {tp_ac}')
plt.legend()
#plt.savefig(f'fig/ac_{name}_gt')
#plt.show()






