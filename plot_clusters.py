import matplotlib.pyplot as plt
import numpy as np

num_clusters = np.load('Predictions/pmi/num_clusters.npy')
num_clusters_ac = np.load('Predictions/ac/num_clusters.npy')
times = np.load('Predictions/pmi/times.npy')
plt.plot(times[:50],num_clusters[:50],label='PMI')
plt.plot(times[:50],num_clusters_ac[:50],label='AC')
plt.xlabel('Markov Time')
plt.ylabel('Number of clusters')
plt.title('Number of Clusters for each Markov Time')
plt.legend()
plt.show()
