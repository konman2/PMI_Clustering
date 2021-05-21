import scipy.io as sio
import numpy as np
from tqdm import tqdm
from sklearn.metrics import normalized_mutual_info_score

taus = 100

filename = f'stability_clustering_tau_1_{taus}_prediction.mat'
predicted = sio.loadmat(filename)['C']

macro_comms = np.load('macro_comms.npy')
micro_comms = np.load('micro_comms.npy')

# Compute the scores.
iterator = tqdm(range(1, 1+taus))
iterator.set_description(f'Computing NMI')
micro_scores = []
macro_scores = []
for tau in iterator:
    macro_score = normalized_mutual_info_score(macro_comms, predicted[:, tau-1])
    micro_score = normalized_mutual_info_score(micro_comms, predicted[:, tau-1])
    macro_scores.append(macro_score)
    micro_scores.append(micro_score)

np.save(f'stability_macro_clustering_tau_1_{taus}.npy', np.asarray(macro_scores))
np.save(f'stability_micro_clustering_tau_1_{taus}.npy', np.asarray(micro_scores))