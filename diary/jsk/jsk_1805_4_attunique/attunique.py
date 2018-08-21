import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import
netfilename = "PKN_22.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

dataname = 'obj_attractor.csv'
conddata = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)

iptidx = []#indicies of the input conditions
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)

num_nodes = len(node_list)
num_edges = len(netdata)

#calculating maximum indegree to set the bounds below
indeglist = {}
for node in node_list:
    indeglist[node] = 0
for link in netdata:
    indeglist[link[2]] += 1
maxindeg = indeglist[max(indeglist, key = indeglist.get)]
wr = 2 * (maxindeg + 2)

paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.float64)
bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.float64)
params = []

for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)
    params.append(np.concatenate((ws[datano], bs[datano])))
params = np.array(params)

tempbs = np.transpose(bs).copy()

bs = np.zeros((len(node_list), totdatano), dtype=np.float64)
i = 0
for idx in range(len(node_list)):
    if idx in iptidx:
        bs[idx] = np.array([99] * totdatano)
    else:
        bs[idx] = tempbs[i]
        i += 1
bs = np.transpose(bs)

netparamslist = []
for w, b in zip(ws, bs):
    netparamslist.append((w, b))

bsdf = pd.DataFrame(bs)
wsdf = pd.DataFrame(ws)


#simulator
signs = np.ones((num_edges,), dtype=np.int32)

ind_rows = np.zeros((num_edges,), dtype=np.int32)
ind_cols = np.zeros((num_edges,), dtype=np.int32)

for idx, edge in enumerate(netdata):
    if edge[1] == 'inhibit':
        signs[idx] = -1
    ind_rows[idx] = node_list.index(edge[0])
    ind_cols[idx] = node_list.index(edge[2])
# end of for
weight_matrix = np.zeros((num_nodes, num_nodes))


def _index_hist(st, hist):
    for i, state in enumerate(hist):
        if (state == st).all():
            return i
    return -1
    
###the Boolean network model
def _func(w, b, inistates):
    weight_matrix[ind_rows, ind_cols] = signs*np.round(w)
    #self.weight_matrix = self.weight_matrix.astype(np.int64)
    #get the attractors
    state_history = [] # Record of state transition
    
    x_t1 = np.array(inistates)
    state_history.append(x_t1)
    while True:#state transition from t to t + 1
        
        x_t2 = np.dot(x_t1, weight_matrix) + b #transitioned state
        
        # Forcing the initial states
        for idx in iptidx:             
            if inistates[idx] == 0:
                x_t2[idx] = 0
            elif inistates[idx] == 1:
                x_t2[idx] = 1
        # end of for

        # Activation function
        x_t2[x_t2 > 0] = 1
        x_t2[x_t2 <= 0] = 0
        
        idx = _index_hist(x_t2, state_history)
        if idx >= 0:
            # point attractor has list length of 1
            attractor = state_history[idx:]
            return attractor
        
        else: # move on
            state_history.append(x_t2.copy())
            x_t1 = x_t2


attractorlist = []
incornodeslist = []
for iws, ibs in netparamslist:
    attractors = []
    incornodes = []
    for condno in range(int(len(conddata)/2)):
        condno += 1
        inistates = conddata.loc['INPUT: ' + str(condno)].tolist()
        for idx in iptidx:#fixing input nodes using extreme basal activity settings
            if inistates[idx] == 0:
                ibs[idx] = -wr
            elif inistates[idx] == 1:
                ibs[idx] = wr
        inistates = np.array(inistates)
        inistates[inistates == 9] = 0
        attractor = _func(iws, ibs, inistates)#calculated states form the Boolean model
        attractors.append(attractor)
        objstates = conddata.loc['OUTPUT: ' + str(condno)].tolist()#objective states
        objidx = [idx for idx, state in enumerate(objstates) if state != 9]
        incornode = []
        for attstate in attractor:
            for idx in objidx:
                if objstates[idx] != attstate[idx]:
                    incornode.append(node_list[idx])
            incornode.append('\t')
        incornodes.append(incornode)
    attractorlist.append(attractors)#attractors for best fit parameters
    incornodeslist.append(incornodes)#incorrect nodes for best fit parameters; none if perfect

attlist = []
for attractors in attractorlist:
    att = np.zeros((3, num_nodes))
    for i, attractor in enumerate(attractors):
        att[i] = np.mean(attractor, axis=0)
    att = att.flatten()
    attlist.append(att)
attlist = np.array(attlist)

models = np.arange(1, len(attlist) + 1)
unimodels = [len(np.unique(np.array(attlist[:i+1]), axis=0)) for i in range(len(attlist))]
plt.figure()
plt.plot(models, unimodels, models, models)









