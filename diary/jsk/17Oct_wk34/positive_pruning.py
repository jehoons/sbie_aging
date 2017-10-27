import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

miraw = pd.read_csv('MI_raw.csv', index_col=0)
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)
#a/f feedback

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
#miraw > 2.4
hicormi24 = miraw > 2.4
n24 = sum(np.array(hicormi24).flatten().astype(int))

print('miraw > 2.4, ' + str(n24))
himi24 = []
for i in range(len(hicormi24)):
    for j in range(len(hicormi24)):
        if hicormi24.iloc[i][j]:
            print(hicormi24.index[i], hicormi24.columns.values[j])
            himi24.append([hicormi24.index[i], hicormi24.columns.values[j], miraw.iloc[i][j]])

#miraw > 2.2
hicormi22 = miraw > 2.2
n22 = sum(np.array(hicormi22).flatten().astype(int))

print('miraw > 2.2, ' + str(n22))
himi22 = []
for i in range(len(hicormi22)):
    for j in range(len(hicormi22)):
        if hicormi22.iloc[i][j]:
            print(hicormi22.index[i], hicormi22.columns.values[j])
            himi22.append([hicormi22.index[i], hicormi22.columns.values[j], miraw.iloc[i][j]])

#miraw > 2.0
hicormi20 = miraw > 2.0
n20 = sum(np.array(hicormi20).flatten().astype(int))

print('miraw > 2.0, ' + str(n20))
himi20 = []
for i in range(len(hicormi20)):
    for j in range(len(hicormi20)):
        if hicormi20.iloc[i][j]:
            print(hicormi20.index[i], hicormi20.columns.values[j])
            himi20.append([hicormi20.index[i], hicormi20.columns.values[j], miraw.iloc[i][j]])

#miraw > 1.8
hicormi18 = miraw > 1.8
n18 = sum(np.array(hicormi18).flatten().astype(int))

print('miraw > 1.8, ' + str(n18))
himi18 = []
for i in range(len(hicormi18)):
    for j in range(len(hicormi18)):
        if hicormi18.iloc[i][j]:
            print(hicormi18.index[i], hicormi18.columns.values[j])
            himi18.append([hicormi18.index[i], hicormi18.columns.values[j], miraw.iloc[i][j]])

#miraw > .1; original
hicormiori = miraw > .1
nori = sum(np.array(hicormiori).flatten().astype(int))

print('pruned with all data, ' + str(nori))
himiori = []
for i in range(len(hicormiori)):
    for j in range(len(hicormiori)):
        if hicormiori.iloc[i][j]:
            print(hicormiori.index[i], hicormiori.columns.values[j])
            himiori.append([hicormiori.index[i], hicormiori.columns.values[j], miraw.iloc[i][j]])

#prune existing PKN using the pruned interaction information
agingnetwork24 = []
for sigmi in himi24:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in agingnetwork.index:
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork24:
                        agingnetwork24.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork24:
                        agingnetwork24.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])

agingnetwork22 = []
for sigmi in himi22:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in agingnetwork.index:
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork22:
                        agingnetwork22.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork22:
                        agingnetwork22.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])

agingnetwork20 = []
for sigmi in himi20:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in agingnetwork.index:
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork20:
                        agingnetwork20.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork20:
                        agingnetwork20.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])

agingnetwork18 = []
for sigmi in himi18:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in agingnetwork.index:
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork18:
                        agingnetwork18.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetwork18:
                        agingnetwork18.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])

agingnetworkori = []
for sigmi in himiori:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in agingnetwork.index:
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetworkori:
                        agingnetworkori.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx]:
                if abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                    if not [agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]] in agingnetworkori:
                        agingnetworkori.append([agingnetwork['Source'][idx], agingnetwork['Regulation Type'][idx], agingnetwork['Target'][idx]])

#export the resulting network files to .sif format
agingnetwork24 = pd.DataFrame(agingnetwork24)
agingnetwork24.to_csv('agingnetwork24p.sif', sep='\t', index=None, header=None)

agingnetwork22 = pd.DataFrame(agingnetwork22)
agingnetwork22.to_csv('agingnetwork22p.sif', sep='\t', index=None, header=None)

agingnetwork20 = pd.DataFrame(agingnetwork20)
agingnetwork20.to_csv('agingnetwork20p.sif', sep='\t', index=None, header=None)

agingnetwork18 = pd.DataFrame(agingnetwork18)
agingnetwork18.to_csv('agingnetwork18p.sif', sep='\t', index=None, header=None)

agingnetworkori = pd.DataFrame(agingnetworkori)
agingnetworkori.to_csv('agingnetworkori.sif', sep='\t', index=None, header=None)















































