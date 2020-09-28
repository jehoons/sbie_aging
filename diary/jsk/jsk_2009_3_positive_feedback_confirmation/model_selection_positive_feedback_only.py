import numpy as np
import pandas as pd
import pickle


#algorithm parameters
n = 5000
positive_feedback = {'PDK1': 1, 'AKT1': 1, 'IKBKB': 1, 'PTEN': 0}
positive_feedback_nodes = ['PDK1', 'AKT1', 'IKBKB', 'PTEN']

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
attractor_df = pd.DataFrame(index=np.arange(len(attractorlist)), columns=node_list)
for _i, _attractor in enumerate(attractorlist):
    attractor_df.loc[_i] = np.mean(_attractor, axis=0)
attractor_avg = attractor_df.mean()

#select positive feedback models
attractor_df_pos_feed = attractor_df[positive_feedback_nodes]
model_idx_pos_feed = []
for _idx, state in attractor_df_pos_feed.iterrows():
    if all([state[pos_feed_node] == positive_feedback[pos_feed_node]for pos_feed_node in positive_feedback_nodes]):
        model_idx_pos_feed.append(_idx)
attractor_df_ = attractor_df.iloc[model_idx_pos_feed]
attractor_df_pos_feed_ = attractor_df_[positive_feedback_nodes]
with open('positive_feedback_model.bin', 'wb') as f:
    pickle.dump(model_idx_pos_feed, f)
























