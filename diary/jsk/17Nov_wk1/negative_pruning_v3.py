import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

###This code prunes aging network using miraw
###If a link present in the agingnetwork has a lower than cutoff mutual information
###It will be removed in the processed network

def negprune(miraw, agingnetwork, abmap, frac):
    #plotting the relationship between mutual information cutoff and the number of interactions
    colist=[]
    nlist=[]
    for cutoff in range(1, 30):
        cutoff = cutoff/10
        hicormi = (miraw < cutoff) & (miraw > 0.000001)
        n = sum(np.array(hicormi).flatten().astype(int))
        colist.append(cutoff)
        nlist.append(n)
    plt.plot(colist, nlist)
    plt.xlabel('cutoff')
    plt.ylabel('MI pairs')
    plt.show()
    
    #selecting the appropriate cutoff; we used the cutoff that allows frac amount of samples to be considered to the tenth
    idx = np.argmin(abs(np.array(nlist) - np.size(np.array(miraw)) * frac / 2))
    cutoff = colist[idx]
    
    #prune selected interactions based on the appropriate cutoffs
    lowcormi = ((miraw < cutoff) & (miraw > .000001))
    himi = []
    for col in lowcormi.columns.values:
        for idx in lowcormi.index:
            if lowcormi[idx][col]:
                himi.append([idx, col, miraw[idx][col]])
    
    #prune existing PKN using the pruned interaction information
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
    
    #return the resulting network file and its cutoff
    return pd.DataFrame(tmpagingnetwork), cutoff


if __name__ == "__main__":
    #importing relevant data
    miraw = pd.read_csv('MI_raw.csv', index_col=0)
    agingnetwork = pd.read_csv('Aging_network.csv')
    abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

    #determine the appropriate agingnetwork
    agnet, cutoff = negprune(miraw, agingnetwork, abmap, .5)

    #export the resulting network files to .sif format
    agnet.to_csv('agingnetwork' + str(int(cutoff*10)) + 'n.sif', sep='\t', index=None, header=None)










