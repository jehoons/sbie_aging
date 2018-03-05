import numpy as np
import pandas as pd

#read netparamsdump.txt along with others
netfilename = "pkn_21_cut0.0.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)

paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

def linkidx(src, tgt, netdata = np.array(netdata)):
    return np.argwhere(np.logical_and(netdata[:,0] == src, netdata[:,2] == tgt))[0][0]

totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.int64)
bs = np.zeros((totdatano, len(node_list)), dtype=np.int64)
for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)

def printbinfo(node):
    print(node)
    print(bs[:, node_list.index(node)])

def printwinfo(src, tgt):
    print(src + ' int w/ ' + tgt)
    print(ws[:, linkidx(src, tgt)])

pnodelist = ['NFKB1', 'SIRT1', 'NFKBIE', 'IL1B', 'IL6', 'TNF']
pedgelist = [['SIRT1', 'NFKB1'], ['NFKBIE', 'NFKB1'], ['NFKB1', 'IL1B'], ['NFKB1', 'IL6'], ['NFKB1', 'TNF']]

for nd in pnodelist:
    printbinfo(nd)
print('\n')
for eg in pedgelist:
    printwinfo(eg[0], eg[1])
