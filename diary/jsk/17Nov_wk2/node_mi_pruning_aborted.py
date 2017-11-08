import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

###This code prunes aging network using miraw
###It selects nodes that are present in both PKN and the dataset
###If an interaction is already defined in PKN, it is maintained,
###and if not, the high Mutual Information links are maintained as it is

#import relavent datasets
miraw = pd.read_csv('MI_raw.csv', index_col=0)
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

#collect nodes that are present in both PKN and dataset
nodeset = []
for ab in abmap.index:
    if abmap['include network'][ab] == 'o' and not abmap['symbol'][ab] in nodeset:
        nodeset.append(abmap['symbol'][ab])

#plotting the relationship between mutual information cutoff and the number of interactions
colist=[]
nlist=[]
for cutoff in range(1, 30):
    cutoff = cutoff/10
    hicormi = miraw > cutoff
    n = sum(np.array(hicormi).flatten().astype(int))
    colist.append(cutoff)
    nlist.append(n)
plt.plot(colist, nlist)
plt.xlabel('cutoff')
plt.ylabel('MI pairs')
plt.show()

#selecting the appropriate cutoff; we used the cutoff that allows frac amount of samples to be considered to the tenth
frac = .1
idx = np.argmin(abs(np.array(nlist) - np.size(np.array(miraw)) * frac / 2))
cutoff = colist[idx]
mimask = miraw > cutoff

#connect the nodes that have above cutoff mutual information according to the miraw data
tempnet = []
for rowab in mimask.index:
    for colab in mimask.columns:
        if mimask[colab][rowab] and abmap['symbol'][colab] in nodeset and abmap['symbol'][rowab] in nodeset:
            if not [abmap['symbol'][colab], 0, abmap['symbol'][rowab]] in tempnet and not [abmap['symbol'][rowab], 0, abmap['symbol'][colab]] in tempnet:
                tempnet.append([abmap['symbol'][colab], 0, abmap['symbol'][rowab]])





