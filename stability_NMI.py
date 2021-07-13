import scipy.io as sio
import numpy as np
from tqdm import tqdm
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score,adjusted_mutual_info_score
import matplotlib.pyplot as plt

def add(filename,desc,taus,times,ax1,ax2,num=70,func=normalized_mutual_info_score,name='airport_ww'):
    predicted = np.load(filename)
    if name == 'airport':
        name= 'airport_ww'
    macro_comms = np.load(f'Graphs/{name}/macro_comms.npy')
    micro_comms = np.load(f'Graphs/{name}/micro_comms.npy')
    iterator = tqdm(range(1, taus+1))
    iterator.set_description(f'Computing NMI')
    micro_scores = []
    macro_scores = []
    for tau in iterator:
        #print(predicted.shape)
        macro_score = func(macro_comms, predicted[:, tau-1])
        micro_score = func(micro_comms, predicted[:, tau-1])
        macro_scores.append(macro_score)
        micro_scores.append(micro_score)
    best_mac = np.argmax(macro_scores)
    best_mic = np.argmax(micro_scores)
    inds = []
    

    ax1.plot(times[:num],macro_scores[:num],label='{} (best is {:.5f} at time {:.5f})'.format(desc,macro_scores[best_mac],times[best_mac]))
    ax2.plot(times[:num],micro_scores[:num],label='{} (best is {:.5f} at time {:.5f})'.format(desc,micro_scores[best_mic],times[best_mic]))


# exit()
def graph(name,taus):

    # filename = f'stability_clustering_tau_1_{taus}_prediction.mat'
    filename = f'Predictions/{name}/pmi/predicted_communities_{taus}.npy'
    filename_ac = f'Predictions/{name}/ac/predicted_communities_{taus}.npy'

    filename_lmepmi = f'Predictions/{name}/lmepmi/predicted_communities_{taus}.npy'
    filename_mac = f'Predictions/{name}/mac/predicted_communities_{taus}.npy'

    times = np.load(f'Predictions/{name}/pmi/times.npy')
    times_oth = np.load(f'Predictions/{name}/lmepmi/times.npy')
    fig1,ax1 = plt.subplots()
    ax1.set_xlabel('Markov time log scale')
    ax1.set_ylabel('NMI')
    ax1.set_title(f'{name} Network Macro Communities')

    fig2,ax2 = plt.subplots()
    ax2.set_xlabel('Markov time log scale')
    ax2.set_ylabel('NMI')
    ax2.set_title(f'{name} Network Micro Communities')
    times = np.log10(times)
    times_oth = np.log10(times_oth)
    func = normalized_mutual_info_score
    add(filename,'PMI',taus,times,ax1,ax2,num=taus,func=func,name=name)
    add(filename_ac,'AC',taus,times,ax1,ax2,num=taus,func=func,name=name)
    add(filename_lmepmi,'lmepmi',taus,times_oth,ax1,ax2,num=taus,func=func,name=name)
    add(filename_mac,'mac',taus,times_oth,ax1,ax2,num=taus,func=func,name=name)
    ax1.legend()
    ax2.legend()
    plt.show()

graph('entsoe',50)

# #predicted = sio.loadmat(filename)['C']
# predicted = np.load(filename)
# predicted_ac = np.load(filename_ac)
# predicted_lmepmi = np.load(filename_lmepmi)

# macro_comms = np.load('Graphs/airport_ww/macro_comms.npy')
# micro_comms = np.load('Graphs/airport_ww/micro_comms.npy')

# # Compute the scores.
# iterator = tqdm(range(1, 1+taus))
# iterator.set_description(f'Computing NMI')
# micro_scores_ac = []
# macro_scores_ac = []
# micro_scores = []
# macro_scores = []
# for tau in iterator:
#     macro_score = normalized_mutual_info_score(macro_comms, predicted[:, tau-1])
#     micro_score = normalized_mutual_info_score(micro_comms, predicted[:, tau-1])
#     macro_scores.append(macro_score)
#     micro_scores.append(micro_score)
#     macro_score_ac = normalized_mutual_info_score(macro_comms, predicted_ac[:, tau-1])
#     micro_score_ac = normalized_mutual_info_score(micro_comms, predicted_ac[:, tau-1])
#     macro_scores_ac.append(macro_score_ac)
#     micro_scores_ac.append(micro_score_ac)
# np.save(f'pmi_macro_clustering_tau_1_{taus}.npy', np.asarray(macro_scores))
# np.save(f'pmi_micro_clustering_tau_1_{taus}.npy', np.asarray(micro_scores))
# np.save(f'ac_macro_clustering_tau_1_{taus}.npy', np.asarray(macro_scores))
# np.save(f'ac_micro_clustering_tau_1_{taus}.npy', np.asarray(micro_scores))
# best_mac = np.argmax(macro_scores)
# best_mac_ac = np.argmax(macro_scores_ac)
# best_mic_ac = np.argmax(micro_scores_ac)
# print("Best Macro Score:",times[best_mac],macro_scores[best_mac])
# best_mic = np.argmax(micro_scores)
# print(best_mic,best_mac)
# print(best_mic_ac,best_mac_ac)
# print("Best Micro Score:",times[best_mic],micro_scores[best_mic])
# plt.figure(1)
# plt.plot(times[:70],macro_scores[:70],label='PMI (best is {:.5f} at time {:.5f})'.format(macro_scores[best_mac],times[best_mac]))
# plt.plot(times[:70],macro_scores_ac[:70],label='AC(best is {:.5f} at time {:.5f})'.format(macro_scores_ac[best_mac_ac],times[best_mac_ac]))
# plt.xlabel('Markov time')
# plt.ylabel('NMI')
# plt.legend()
# plt.title('Airport Network Macro Communities')
# plt.figure(2)
# plt.plot(times[:70],micro_scores[:70],label='PMI (best is {:.5f} at time {:.5f})'.format(micro_scores[best_mic],times[best_mic]))
# plt.plot(times[:70],micro_scores_ac[:70],label='AC (best is {:.5f} at time {:.5f})'.format(micro_scores_ac[best_mic_ac],times[best_mic_ac]))
# plt.xlabel('Markov time')
# plt.ylabel('NMI')
# plt.legend()
# plt.title('Airport Network Micro Communities')
# plt.show()
