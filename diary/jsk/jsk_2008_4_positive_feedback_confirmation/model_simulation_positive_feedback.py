import numpy as np
import pandas as pd
import itertools as it


#algorithm parameters
n = 5000

#import network, network parameter and input-output conditions
netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

conddata = pd.read_csv('obj_attractor_control_ss.csv', index_col = 0)
node_list = list(conddata.columns)

#generate output conditions
objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
objboolmask = objstates != 9
objphen = objstates[objboolmask]

iptidx = []#indicies of the input conditions
iptval = []#input value for code revision
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)
        if state == 1:
            iptval.append(99)
        else:
            iptval.append(-99)

optidx =[]#indicies of the output conditions
for idx, state in enumerate(list(conddata.iloc[1])):
    if state != 9:
        optidx.append(idx)

num_nodes = len(node_list)
num_edges = len(netdata)

netparamslist = []

paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.float64)
bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.float64)

for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)

tempbs = np.transpose(bs).copy()

bs = np.zeros((len(node_list), totdatano), dtype=np.float64)
i = 0
for idx in range(len(node_list)):
    if idx in iptidx:
        bs[idx] = np.array([iptval[iptidx.index(idx)]] * totdatano)
    else:
        bs[idx] = tempbs[i]
        i += 1
bs = np.transpose(bs)

for w, b in zip(ws, bs):
    netparamslist.append((w, b))
netparamslist = netparamslist[:n]


#run Boolean network simulation
def _index_hist(st, hist):
    for i, state in enumerate(hist):
        if (state == st).all():
            return i
    return -1

attractorlist = []
statehistorylist = []
for ws, bs in netparamslist:
    _netdata = netdata.copy()
    _node_list = node_list.copy()
    
    signs = np.ones((num_edges,), dtype=np.int32)
    ind_rows = np.zeros((num_edges,), dtype=np.int32)
    ind_cols = np.zeros((num_edges,), dtype=np.int32)
    
    for idx, edge in enumerate(_netdata):
        if edge[1] == 'inhibit':
            signs[idx] = -1
        ind_rows[idx] = _node_list.index(edge[0])
        ind_cols[idx] = _node_list.index(edge[2])
    # end of for
    
    ###the Boolean network model
    def _func(b, inistates):
        weight_matrix = np.zeros((len(_node_list), len(_node_list)))
        weight_matrix[ind_rows, ind_cols] = signs*np.round(ws)
        
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
                return attractor, state_history
            
            else: # move on
                state_history.append(x_t2.copy())
                x_t1 = x_t2
    
    #set initial conditions
    inistates = conddata.loc['INPUT'].tolist()
    for idx in iptidx:#fixing input nodes using extreme basal activity settings
        if inistates[idx] == 0:
            bs[idx] = -99
        elif inistates[idx] == 1:
            bs[idx] = 99
    inistates = np.array(inistates)
    inistates[inistates == 9] = 0
    
    attractor, statehistory = _func(bs, inistates)
    attractorlist.append(attractor)
    statehistorylist.append(statehistory)


#analyze resulting attractor information
attractor_df = pd.DataFrame(index=np.arange(n), columns=node_list)
for _i, _attractor in enumerate(attractorlist):
    attractor_df.loc[_i] = np.mean(_attractor, axis=0)
attractor_avg = attractor_df.mean()


#attractor state after senesence input is removed
attractorlist_new = []
statehistorylist_new = []
for (ws, bs), attractor_prev in zip(netparamslist, attractorlist):
    _netdata = netdata.copy()
    _node_list = node_list.copy()
    
    signs = np.ones((num_edges,), dtype=np.int32)
    ind_rows = np.zeros((num_edges,), dtype=np.int32)
    ind_cols = np.zeros((num_edges,), dtype=np.int32)
    
    for idx, edge in enumerate(_netdata):
        if edge[1] == 'inhibit':
            signs[idx] = -1
        ind_rows[idx] = _node_list.index(edge[0])
        ind_cols[idx] = _node_list.index(edge[2])
    # end of for
    
    ###the Boolean network model
    def _func(b, inistates):
        weight_matrix = np.zeros((len(_node_list), len(_node_list)))
        weight_matrix[ind_rows, ind_cols] = signs*np.round(ws)
        
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
                return attractor, state_history
            
            else: # move on
                state_history.append(x_t2.copy())
                x_t1 = x_t2
    
    #set initial conditions to previous senescence attractor with input signal turned off
    inistates = np.mean(attractor_prev, axis=0)
    for idx in iptidx:#fixing input nodes using extreme basal activity settings
        inistates[idx] = 0
        bs[idx] = -99
    
    attractor, statehistory = _func(bs, inistates)
    attractorlist_new.append(attractor)
    statehistorylist_new.append(statehistory)


#analyze resulting attractor information
attractor_df_new = pd.DataFrame(index=np.arange(n), columns=node_list)
for _i, _attractor in enumerate(attractorlist_new):
    attractor_df_new.loc[_i] = np.mean(_attractor, axis=0)
attractor_avg_new = attractor_df_new.mean()


#analyze PDK1 positive feedback
#positive_feedback = {'PDK1': 1, 'AKT': 1, 'IKBKB': 1, 'PTEN': 0}
positive_feedback_nodes = ['PDK1', 'AKT1', 'IKBKB', 'PTEN']
positive_feedback_bf = attractor_avg.loc[positive_feedback_nodes]
positive_feedback_af = attractor_avg_new.loc[positive_feedback_nodes]






































