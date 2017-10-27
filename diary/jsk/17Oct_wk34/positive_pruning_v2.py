import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#importing relevant data
miraw = pd.read_csv('MI_raw.csv', index_col=0)
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

###This code prunes aging network using miraw
###It finds the intersection of links present in the PKN and high mutual information links from miraw
###Only the links that are well-represented with high mutual information will be present in the processed network

#plotting the relationship between mutual information cutoff and the number of interactions
colist=[]
nlist=[]
for cutoff in range(12, 27):
    cutoff = cutoff/10
    hicormi = miraw > cutoff
    n = sum(np.array(hicormi).flatten().astype(int))
    colist.append(cutoff)
    nlist.append(n)
plt.plot(colist, nlist)
plt.xlabel('cutoff')
plt.ylabel('MI pairs')

#prune selected interactions based on the appropriate cutoffs
cutofflist = [2.4, 2.2, 2.0, 1.8, .1]#list of miraw cutoffs to test; only links above cutoff mi are selected; .1 is for original PKN
hicormilist = []
for cutoff in cutofflist:
    hicormilist.append(miraw > cutoff)
himilist = []
for hicormi in hicormilist:
    himi = []
    for col in hicormi.columns.values:
        for idx in hicormi.index:
            if hicormi[idx][col]:
                himi.append([idx, col, miraw[idx][col]])
    himilist.append(himi)

#prune existing PKN using the pruned interaction information
agingnetworklist = []
for himi in himilist:
    tmpagingnetwork = []
    for sigmi in himi:
        if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
            for idx in agingnetwork.index:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in tmpagingnetwork:
                        tmpagingnetwork.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
                elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in tmpagingnetwork:
                        tmpagingnetwork.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
    agingnetworklist.append(tmpagingnetwork)

#export the resulting network files to .sif format
for cutoffidx, agingnetwork in enumerate(agingnetworklist):
    agingnetwork = pd.DataFrame(agingnetwork)
    agingnetwork.to_csv('agingnetwork' + str(int(cutofflist[cutoffidx]*10)) + 'p.sif', sep='\t', index=None, header=None)

















