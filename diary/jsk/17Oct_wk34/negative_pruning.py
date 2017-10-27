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
    hicormi = (miraw < cutoff) & (miraw > 0.000001)
    n = sum(np.array(hicormi).flatten().astype(int))
    colist.append(cutoff)
    nlist.append(n)
plt.plot(colist, nlist)
plt.xlabel('cutoff')
plt.ylabel('MI pairs')

#prune selected interactions based on the appropriate cutoffs
#miraw < 1.2
hicormi12 = (miraw < 1.2) & (miraw > 0.000001)
n12 = sum(np.array(hicormi12).flatten().astype(int))

print('miraw < 1.2, ' + str(n12))
himi12 = []
for i in range(len(hicormi12)):
    for j in range(len(hicormi12)):
        if hicormi12.iloc[i][j]:
            print(hicormi12.index[i], hicormi12.columns.values[j])
            himi12.append([hicormi12.index[i], hicormi12.columns.values[j], miraw.iloc[i][j]])

#miraw < 1.4
hicormi14 = (miraw < 1.4) & (miraw > 0.000001)
n14 = sum(np.array(hicormi14).flatten().astype(int))

print('miraw < 1.4, ' + str(n14))
himi14 = []
for i in range(len(hicormi14)):
    for j in range(len(hicormi14)):
        if hicormi14.iloc[i][j]:
            print(hicormi14.index[i], hicormi14.columns.values[j])
            himi14.append([hicormi14.index[i], hicormi14.columns.values[j], miraw.iloc[i][j]])

#miraw < 1.6
hicormi16 = (miraw < 1.6) & (miraw > 0.000001)
n16 = sum(np.array(hicormi16).flatten().astype(int))

print('miraw < 1.6, ' + str(n16))
himi16 = []
for i in range(len(hicormi16)):
    for j in range(len(hicormi16)):
        if hicormi16.iloc[i][j]:
            print(hicormi16.index[i], hicormi16.columns.values[j])
            himi16.append([hicormi16.index[i], hicormi16.columns.values[j], miraw.iloc[i][j]])

#miraw < 1.8
hicormi18 = (miraw < 1.8) & (miraw > 0.000001)
n18 = sum(np.array(hicormi18).flatten().astype(int))

print('miraw < 1.8, ' + str(n18))
himi18 = []
for i in range(len(hicormi18)):
    for j in range(len(hicormi18)):
        if hicormi18.iloc[i][j]:
            print(hicormi18.index[i], hicormi18.columns.values[j])
            himi18.append([hicormi18.index[i], hicormi18.columns.values[j], miraw.iloc[i][j]])

#miraw < 2.0
hicormi20 = (miraw < 2.0) & (miraw > 0.000001)
n20 = sum(np.array(hicormi20).flatten().astype(int))

print('miraw < 2.0, ' + str(n20))
himi20 = []
for i in range(len(hicormi20)):
    for j in range(len(hicormi20)):
        if hicormi20.iloc[i][j]:
            print(hicormi20.index[i], hicormi20.columns.values[j])
            himi20.append([hicormi20.index[i], hicormi20.columns.values[j], miraw.iloc[i][j]])

#miraw < 2.2
hicormi22 = (miraw < 2.2) & (miraw > 0.000001)
n22 = sum(np.array(hicormi22).flatten().astype(int))

print('miraw < 2.2, ' + str(n22))
himi22 = []
for i in range(len(hicormi22)):
    for j in range(len(hicormi22)):
        if hicormi22.iloc[i][j]:
            print(hicormi22.index[i], hicormi22.columns.values[j])
            himi22.append([hicormi22.index[i], hicormi22.columns.values[j], miraw.iloc[i][j]])

#prune existing PKN using the pruned interaction information
agingnetwork12 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork12.copy()
for sigmi in himi12:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork12:
                    agingnetwork12.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork12:
                    agingnetwork12.remove(temp[idx])

agingnetwork14 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork14.copy()
for sigmi in himi14:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork14:
                    agingnetwork14.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork14:
                    agingnetwork14.remove(temp[idx])

agingnetwork16 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork16.copy()
for sigmi in himi16:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork16:
                    agingnetwork16.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork16:
                    agingnetwork16.remove(temp[idx])

agingnetwork18 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork18.copy()
for sigmi in himi18:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork18:
                    agingnetwork18.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork18:
                    agingnetwork18.remove(temp[idx])

agingnetwork20 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork20.copy()
for sigmi in himi20:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork20:
                    agingnetwork20.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork20:
                    agingnetwork20.remove(temp[idx])

agingnetwork22 = np.array([agingnetwork['Source'].values.tolist(),agingnetwork['Regulation Type'].values.tolist(),agingnetwork['Target'].values.tolist()]).T.tolist()
temp = agingnetwork22.copy()
for sigmi in himi22:
    if abmap['include network'][sigmi[0]] == 'o' and abmap['include network'][sigmi[1]] == 'o':
        for idx in range(len(temp)):
            if abmap['symbol'][sigmi[0]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[1]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork22:
                    agingnetwork22.remove(temp[idx])
            elif abmap['symbol'][sigmi[1]] == agingnetwork['Source'][idx] and abmap['symbol'][sigmi[0]] == agingnetwork['Target'][idx]:
                if temp[idx] in agingnetwork22:
                    agingnetwork22.remove(temp[idx])

#export the resulting network files to .sif format
agingnetwork12 = pd.DataFrame(agingnetwork12)
agingnetwork12.to_csv('agingnetwork12n.sif', sep='\t', index=None, header=None)

agingnetwork14 = pd.DataFrame(agingnetwork14)
agingnetwork14.to_csv('agingnetwork14n.sif', sep='\t', index=None, header=None)

agingnetwork16 = pd.DataFrame(agingnetwork16)
agingnetwork16.to_csv('agingnetwork16n.sif', sep='\t', index=None, header=None)

agingnetwork18 = pd.DataFrame(agingnetwork18)
agingnetwork18.to_csv('agingnetwork18n.sif', sep='\t', index=None, header=None)

agingnetwork20 = pd.DataFrame(agingnetwork20)
agingnetwork20.to_csv('agingnetwork20n.sif', sep='\t', index=None, header=None)

agingnetwork22 = pd.DataFrame(agingnetwork22)
agingnetwork22.to_csv('agingnetwork22n.sif', sep='\t', index=None, header=None)















































