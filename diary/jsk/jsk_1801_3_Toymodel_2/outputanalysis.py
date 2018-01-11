import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#importing data to analyze
netno = 8
linklist = []
with open('PKN_20_' + str(netno) + '.sif', 'r') as netfile:
    for link in netfile:
        linklist.append(link.strip().split('\t'))
conddata = pd.read_csv('obj_attractor.csv', index_col = 0)
nodelist = []
for idx, state in enumerate(list(conddata.iloc[0])):
    if state == 9:
        nodelist.append(conddata.columns[idx])
data = []
with open('output_fittedparams_PKN_20_' + str(netno) + '.txt', 'r') as f:
    while True:
        line = f.readline()
        data.append(line)
        if not line: break

#reading the output data file
fitlist = []
wslist = []
bslist = []
for idx, line in enumerate(data):
    if line == 'Weights:\n':
        wslist.append(data[idx + 1].strip().split('\t'))
    elif line == 'Basal Activities:\n':
        bslist.append(data[idx + 1].strip().split('\t'))
    else:
        try:
            stuff = line.strip().split(' ')
            if stuff[1] == 'Best':
                fitlist.append(float(stuff[3]))
        except IndexError:
            continue
wslist = np.array(wslist, dtype = np.int32)
bslist = np.array(bslist, dtype = np.int32)

#only select results with fitness 0
pwslist = []
pbslist = []
for fit, ws, bs in zip(fitlist, wslist, bslist):
    if fit > 0:
        continue
    else:
        pwslist.append(ws)
        pbslist.append(bs)

pwslistt = np.array(pwslist).T
pbslistt = np.array(pbslist).T


#draw histogram of each weight
ncol = 4#number of histograms to fit in a row
wbin = range(np.amin(pwslistt), np.amax(pwslistt) + 2)
fig, subp = plt.subplots(int(len(pwslistt) / ncol + int(len(pwslistt) % ncol > 0)), ncol, sharex = True)
for n, (pws, link) in enumerate(zip(pwslistt, linklist)):
    subp[int(n / ncol), n % ncol].hist(pws, bins = wbin, density = True)
    subp[int(n / ncol), n % ncol].set_title(link)

#draw histogram of each weight
ncol = 3#number of histograms to fit in a row
bbin = range(np.amin(pbslistt), np.amax(pbslistt) + 2)
fig, subp = plt.subplots(int(len(pbslistt) / ncol + int(len(pbslistt) % ncol > 0)), ncol, sharex = True)
for n, (pbs, node) in enumerate(zip(pbslistt, nodelist)):
    subp[int(n / ncol), n % ncol].hist(pbs, bins = bbin, density = True)
    subp[int(n / ncol), n % ncol].set_title(node)

#draw pairwise scatterplot of basal activities and weights
fig, subp = plt.subplots(len(nodelist), len(nodelist))
bbin = range(np.amin(pbslistt), np.amax(pbslistt) + 2)
wmin = np.amin(pbslistt)
wmax = np.amax(pbslistt)
for row in range(len(nodelist)):
    for col in range(len(nodelist)):
        if row == col:
            subp[row, col].hist(pbslistt[row], bins = bbin, density = True)
        else:
            subp[row, col].scatter(pbslistt[col], pbslistt[row])#misimplemented
            subp[row, col].set_xlim(wmin, wmax)
            subp[row, col].set_ylim(wmin, wmax)



