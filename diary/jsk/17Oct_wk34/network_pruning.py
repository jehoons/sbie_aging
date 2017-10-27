import pandas as pd
import numpy as np

miraw = pd.read_csv('MI_raw.csv', index_col=0)
mis = pd.read_csv('MI_s.csv', index_col=0)
miq = pd.read_csv('MI_q.csv', index_col=0)
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

#b/f feedback
hicormi26 = miraw > 2.6
n26=sum(np.array(hicormi26).flatten().astype(int))#20

hicormi255 = miraw > 2.55
n255=sum(np.array(hicormi255).flatten().astype(int))#109

hicormis138 = mis > 1.38
ns138=sum(np.array(hicormis138).flatten().astype(int))#5151 - not enough data; 0 at 1.39

hicormiq138 = miq > 1.38
nq138=sum(np.array(hicormiq138).flatten().astype(int))#22155 - not enough data; 0 at 1.39

print('miraw > 2.6, ' + str(n26))
himi26 = []
for i in range(len(hicormi26)):
    for j in range(len(hicormi26)):
        if hicormi26.iloc[i][j]:
            print(hicormi26.index[i], hicormi26.columns.values[j], miraw.iloc[i][j])
            himi26.append([hicormi26.index[i], hicormi26.columns.values[j], miraw.iloc[i][j]])

print('miraw > 2.55, ' + str(n255))
himi255 = []
for i in range(len(hicormi255)):
    for j in range(len(hicormi255)):
        if hicormi255.iloc[i][j]:
            print(hicormi255.index[i], hicormi255.columns.values[j], miraw.iloc[i][j])
            himi255.append([hicormi255.index[i], hicormi255.columns.values[j], miraw.iloc[i][j]])






