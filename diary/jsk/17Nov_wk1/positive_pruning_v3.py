import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

###This code prunes aging network using miraw
###It finds the intersection of links present in the PKN and high mutual information links from miraw
###Only the links that are well-represented with high mutual information will be present in the processed network

def posprune(miraw, agingnetwork, abmap):
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
    
    #selecting the appropriate cutoff; we used the cutoff closest to the median number of interactions to the tenth
    idx = np.argmin(abs(np.array(nlist) - np.size(np.array(miraw)) / 4))
    cutoff = colist[idx]
    
    #prune selected interactions based on the appropriate cutoff
    hicormi = miraw > cutoff
    himi = []
    for col in hicormi.columns.values:
        for idx in hicormi.index:
            if hicormi[idx][col]:
                himi.append([idx, col, miraw[idx][col]])
    
    #prune existing PKN using the pruned interaction information
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
    
    #return the resulting network file and its cutoff
    return pd.DataFrame(tmpagingnetwork), cutoff


if __name__ == "__main__":
    #importing relevant data
    miraw = pd.read_csv('MI_raw.csv', index_col=0)
    agingnetwork = pd.read_csv('Aging_network.csv')
    abmap = pd.read_csv('antibody_mapping.csv', index_col=0)
    
    #determine the appropriate agingnetwork
    agnet, cutoff = posprune(miraw, agingnetwork, abmap)
    
    #export the resulting network files to .sif format
    agnet.to_csv('agingnetwork' + str(int(cutoff*10)) + 'p.sif', sep='\t', index=None, header=None)












