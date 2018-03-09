import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from optweightandbasal import OptWeightandBasal

#read netparamsdump.txt along with others
netfilename = "PKN_21.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)
iptidx = []#indicies of the input conditions
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)

paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

#calculating maximum indegree to set the bounds below
indeglist = {}
for node in node_list:
    indeglist[node] = 0
for link in netdata:
    indeglist[link[2]] += 1
maxindeg = indeglist[max(indeglist, key = indeglist.get)]
wr = 2 * (maxindeg + 2)

def linkidx(src, tgt, netdata = np.array(netdata)):
    return np.argwhere(np.logical_and(netdata[:,0] == src, netdata[:,2] == tgt))[0][0]

totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.int64)
bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.int64)

for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)

ws = np.transpose(ws)
tempbs = np.transpose(bs)

bs = np.zeros((len(node_list), totdatano), dtype=np.int64)

_, datalen = np.shape(ws)
i = 0
for idx in range(len(node_list)):
    if idx in iptidx:
        bs[idx] = np.zeros(datalen)
    else:
        bs[idx] = tempbs[i]
        i += 1

mw = np.mean(ws, axis=1)
mb = np.mean(bs, axis=1)

prob = OptWeightandBasal()#dim = edgeno + nodeno; weightno + basalno

attractors = []
incornodes = []
fitness = 0.0
for condno in range(prob.num_conds):
    condno += 1
    inistates = prob.conddata.loc['INPUT: ' + str(condno)].tolist()

    inistates = np.array(inistates)
    inistates[inistates == 9] = 0
    attractor = prob._func(mw, mb, inistates)#calculated states form the Boolean model
    attractors.append(attractor)
    objstates = prob.conddata.loc['OUTPUT: ' + str(condno)].tolist()#objective states
    objidx = [idx for idx, state in enumerate(objstates) if state != 9]
    incornode = []
    for idx in objidx:
        if objstates[idx] != attractor[0][idx]:
            incornode.append(prob.node_list[idx])
    incornodes.append(incornode)

print('Total number of models: %d' % datalen)
print('Weights:')
print(mw)
print('Basal Activities:')
print(mb)
print('Attractors:')
for attractor in attractors:
    print(attractor)
print('Incorrect nodes:')
print(incornodes)

with open('the_ensemble_model.txt', 'wt') as f:#documenting results
    f.write('Format: tab separated\nWeights\nBasal Activities\n\n')
    for w in mw:
        f.write('%f\t' % w)
    f.write('\n')
    for b in mb:
        f.write('%f\t' % b)
    f.write('\n\n')

'''
n = len(node_list)
#marking
for idx in range(1, n + 1):
    plt.subplot(n, n, idx)
    plt.title(node_list[idx-1])

#plot basal activities
for idx, b in enumerate(bs):
    idx += 1
    plt.subplot(n, n, idx + (idx - 1) * n)
    plt.hist(b, range=(-wr, wr))

#plot weights
for w, link in zip(ws, netdata):
    src = node_list.index(link[0]) + 1
    tgt = node_list.index(link[2]) + 1
    plt.subplot(n, n, src + (tgt - 1) * n)
    plt.hist(w, range=(-wr, wr))
'''
