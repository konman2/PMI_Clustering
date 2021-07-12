import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import geopandas as gpd
import geoplot as gplt
from shapely.geometry import Point
from matplotlib.cm import get_cmap
from scipy.stats import mode

data = pd.read_csv('Graphs/airport_ww/node_info.csv')

#print(data)
lat = data['Latitude'].to_numpy()
long = data['Longitude'].to_numpy()
world = gpd.read_file(
    gpd.datasets.get_path('naturalearth_lowres')
)
points = data.apply(
    lambda srs: Point(float(srs['Longitude']), float(srs['Latitude'])),
    axis='columns'
)
data_geocoded =  gpd.GeoDataFrame(data, geometry=points)
data_geocoded_pmi =  gpd.GeoDataFrame(data, geometry=points)
data_geocoded_ac =  gpd.GeoDataFrame(data, geometry=points)
print(data_geocoded)
#print(data)
#49, 58
#39,52
taus = 100
filename = f'Predictions/pmi/predicted_communities_{taus}.npy'
filename_ac = f'Predictions/ac/predicted_communities_{taus}.npy'
predictions = np.load(filename)
predictions_ac = np.load(filename_ac)
clusters_micro = predictions[:,49]
clusters_macro = predictions[:,58]
clusters_micro_ac = predictions_ac[:,39]
clusters_macro_ac = predictions_ac[:,52]
macro_comms = np.load('Graphs/airport_ww/macro_comms.npy')
micro_comms = np.load('Graphs/airport_ww/micro_comms.npy')
m = len(np.unique(clusters_macro))
print(m)
#print(list(micro_comms))
#colors_micro = ['o' if clusters_micro[i] != micro_comms[i] else 'w' for i in range(len(clusters_micro))] 
#colors_micro_ac = ['b' if clusters_micro_ac[i] != micro_comms[i] else 'w' for i in range(len(clusters_micro))] 
#print(colors_micro_ac)
inds = []
inds_ac = []
for i in np.unique(micro_comms):
    cluster = np.nonzero(micro_comms==i)
    pred = clusters_micro[cluster]
    pred2 = clusters_micro_ac[cluster]
    m = mode(pred)[0][0]
    wrong = np.nonzero(pred!=m)
    # print(cluster[0])
    # print(wrong[0])
    if len(wrong[0]) > 0:
        inds += list(cluster[0][wrong[0]])

    m2 = mode(pred2)[0][0]
    wrong2 = np.nonzero(pred2!=m2)
    # print(cluster[0])
    # print(wrong[0])
    if len(wrong2[0]) > 0:
        inds_ac += list(cluster[0][wrong2[0]])
print(len(inds),len(inds_ac))



    

    
# for c,i in enumerate(clusters_micro):
#     if i != micro_comms[c]:
#         inds.append(c)
# for c,i in enumerate(clusters_micro_ac):
#     if i != micro_comms[c]:
#         inds_ac.append(c)


data_geocoded_pmi = data_geocoded.iloc[inds]
data_geocoded_ac = data_geocoded.iloc[inds_ac]
#print(len(inds),len(inds_ac))
# print(cities)
# print(world.crs)
# print(list(clusters_micro))
# print(list(micro_comms))
# print(clusters_micro_ac)
ax = gplt.polyplot(world,figsize=(8,5))
ax2 = gplt.polyplot(world,figsize=(8,5))

gplt.pointplot(data_geocoded_pmi, ax=ax,s=1,hue=[1 for i in inds],cmap='Set1',label='pmi mistakes')

gplt.pointplot(data_geocoded_ac, ax=ax,s=1,hue=[2 for i in inds_ac],cmap='tab10',alpha=0.5,label='ac mistakes')

#gplt.pointplot(data_geocoded, ax=ax,s=1,hue=colors_micro_ac)

#gplt.pointplot([p])
plt.legend()
plt.show()
print(len(lat))