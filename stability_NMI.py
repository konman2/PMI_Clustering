import scipy.io as sio
import numpy as np
from tqdm import tqdm
from sklearn.metrics import normalized_mutual_info_score,adjusted_rand_score,adjusted_mutual_info_score
import matplotlib.pyplot as plt

def add(filename,desc,taus,times,ax1,ax2,num=70,func=normalized_mutual_info_score,name='airport_ww'):
    predicted = np.load(filename)
    print(predicted)
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

    # filename_stoch = f'Predictions_stoch/{name}/pmi/predicted_communities_{taus}.npy'
    # filename_ac_stoch = f'Predictions_stoch/{name}/ac/predicted_communities_{taus}.npy'
    # filename_lmepmi_stoch = f'Predictions_stoch/{name}/lmepmi/predicted_communities_{taus}.npy'
    # filename_mac_stoch = f'Predictions_stoch/{name}/mac/predicted_communities_{taus}.npy'
    
    filename_nolog = f'Predictions/{name}/nolog/predicted_communities_{taus}.npy'

    times = np.load(f'Predictions/{name}/pmi/times.npy')
    times_oth = np.load(f'Predictions/{name}/lmepmi/times.npy')
    #times2 = np.load(f'Predictions_stoch/{name}/pmi/times.npy')
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
    #times2 = np.log10(times2)
    func = normalized_mutual_info_score
    add(filename,'PMI',taus,times,ax1,ax2,num=taus,func=func,name=name)
    add(filename_ac,'AC',taus,times,ax1,ax2,num=taus,func=func,name=name)
    # add(filename_lmepmi,'lmepmi',taus,times_oth,ax1,ax2,num=taus,func=func,name=name)
    # add(filename_mac,'mac',taus,times_oth,ax1,ax2,num=taus,func=func,name=name)
    #add(filename_nolog,'direct difference',taus,times,ax1,ax2,num=taus,func=func,name=name)
    # add(filename_stoch,'PMI Stoch',100,times2,ax1,ax2,num=taus,func=func,name=name)
    # add(filename_ac_stoch,'AC Stoch',100,times2,ax1,ax2,num=taus,func=func,name=name)
    
    ax1.legend()
    ax2.legend()
    #print(times2,len(times2))
    plt.show()
    # fig1.savefig(f'../Week 3/{name}_mac.png',dpi=100)
    # fig2.savefig(f'../Week 3/{name}_mic.png',dpi=100)

graph('wiki-fields',50)
