import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#importing relevant data
miraw = pd.read_csv('MI_raw.csv', index_col=0)
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

###This code prunes aging network using miraw
###If a link present in the agingnetwork has a lower than cutoff mutual information
###It will be removed in the processed network

#plotting the relationship between mutual information cutoff and the number of interactions
colist=[]
nlist=[]
for cutoff in range(12, 27):
    cutoff = cutoff/10
    hicormi = (miraw < cutoff) & (miraw > 0.000001)
    n = sum(np.array(hicormi).flatten().astype(int))
    colist.append(cutoff)
    nlist.append(n)
plt.plot(colist, nlist)
plt.xlabel('cutoff')
plt.ylabel('MI pairs')

#prune selected interactions based on the appropriate cutoffs
cutofflist = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2]#list of miraw cutoffs to test; mi less than cutoff are selected for removal
lowcormilist = []
for cutoff in cutofflist:
    lowcormilist.append((miraw < cutoff) & (miraw > .000001))
himilist = []
for lowcormi in lowcormilist:
    himi = []
    for col in lowcormi.columns.values:
        for idx in lowcormi.index:
            if lowcormi[idx][col]:
                himi.append([idx, col, miraw[idx][col]])
    himilist.append(himi)

#prune existing PKN using the pruned interaction information
agingnetworklist = []
for himi in himilist:
    tmpagingnetwork = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
    temp = tmpagingnetwork.copy()
    for sigmi in himi:
        if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
            for idx in range(len(temp)):
                if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if temp[idx] in tmpagingnetwork:
                        tmpagingnetwork.remove(temp[idx])
                elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if temp[idx] in tmpagingnetwork:
                        tmpagingnetwork.remove(temp[idx])
    agingnetworklist.append(tmpagingnetwork)

#export the resulting network files to .sif format
for cutoffidx, agingnetwork in enumerate(agingnetworklist):
    agingnetwork = pd.DataFrame(agingnetwork)
    agingnetwork.to_csv('agingnetwork' + str(int(cutofflist[cutoffidx]*10)) + 'n.sif', sep='\t', index=None, header=None)
















