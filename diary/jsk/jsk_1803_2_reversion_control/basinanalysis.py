import numpy as np
import pandas as pd

#number of random initial state generation
iteration = 1000

#import
netfilename = "PKN_22.sif"
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

num_nodes = len(node_list)
num_edges = len(netdata)

ensemblemodel = "the_ensemble_model.txt"
modeldata = []
modelfile = open(ensemblemodel, 'r')
for params in modelfile:
    modeldata.append(params.strip().split('\t'))

ws = np.array(modeldata[4], dtype=np.float64)
bs = np.array(modeldata[5], dtype=np.float64)

signs = np.ones((num_edges,), dtype=np.int32)
ind_rows = np.zeros((num_edges,), dtype=np.int32)
ind_cols = np.zeros((num_edges,), dtype=np.int32)

for idx, edge in enumerate(netdata):
    if edge[1] == 'inhibit':
        signs[idx] = -1
    ind_rows[idx] = node_list.index(edge[0])
    ind_cols[idx] = node_list.index(edge[2])
# end of for

def _index_hist(st, hist):
    for i, state in enumerate(hist):
        if (state == st).all():
            return i
    return -1

###the Boolean network model
def _func(b, inistates):
    weight_matrix = np.zeros((num_nodes, num_nodes))
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
            return attractor
        
        else: # move on
            state_history.append(x_t2.copy())
            x_t1 = x_t2

num_conds = int(len(conddata.index) / 2)
fitlist = []
for itno in range(iteration):
    fitness = 0.0
    for condno in range(num_conds):#for all of the conditions presented calculate the cost function and add all
        condno += 1
        inistates = conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
        inistates = np.array(inistates)
        inistates[inistates == 9] = np.random.randint(0, 2, len(inistates[inistates == 9]))
        predstateslist = _func(bs, inistates) #calculated states form the Boolean model
        objstates = np.array(conddata.loc['OUTPUT: ' + str(condno)].tolist())#objective states
        objboolmask = objstates != 9
        objphen = objstates[objboolmask]
        tmpfitlist = np.zeros((len(predstateslist),))
        for i, predstates in enumerate(predstateslist):
            predstates = np.array(predstates)
            predphen = predstates[objboolmask]
            tmpfit = np.sum((objphen - predphen)**2)
            tmpfitlist[i] = tmpfit
        fitness += np.mean(tmpfitlist)
    fitlist.append(fitness)




