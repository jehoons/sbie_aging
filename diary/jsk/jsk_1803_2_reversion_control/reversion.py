import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt

#algorithm parameters
mode = 0#basically a Boolean value; if true, algorithm considers a single ensemble model, and if false, algorithm considers multiple fitted models
sign = -1#direction of control; -1 or 1
maxdrugtgt = 3#maximum number of drug used

#import
netfilename = "PKN_22.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

dataname = 'obj_attractor_control_qs.csv'
conddata = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)

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

exidx = np.concatenate((iptidx, optidx))
conlist = np.array(node_list.copy())
conlist = np.delete(conlist, exidx)#list of controllable nodes

num_nodes = len(node_list)
num_edges = len(netdata)

netparamslist = []
if mode:
    ###########
    ensemblemodel = "the_ensemble_model.txt"
    modeldata = []
    modelfile = open(ensemblemodel, 'r')
    for params in modelfile:
        modeldata.append(params.strip().split('\t'))
    
    ws = np.array(modeldata[4], dtype=np.float64)
    bs = np.array(modeldata[5], dtype=np.float64)
    
    netparamslist.append((ws, bs))
    ##########
else:
    ##########
    paramsdumpname = "netparamsdump.txt"
    paramsdata = []
    paramsfile = open(paramsdumpname, 'r')
    for params in paramsfile:
        paramsdata.append(params.strip().split('\t'))
    
    totdatano = int((len(paramsdata) - 4) / 3)
    ws = np.zeros((totdatano, len(netdata)), dtype=np.int64)
    bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.int64)
    
    for datano in range(totdatano):
        ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
        bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)

    tempbs = np.transpose(bs).copy()
    
    bs = np.zeros((len(node_list), totdatano), dtype=np.int64)
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
    ##########

tgtlists = []
for ws, bs in netparamslist:
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
    
    inistates = conddata.loc['INPUT'].tolist()
    for idx in iptidx:#fixing input nodes using extreme basal activity settings
        if inistates[idx] == 0:
            bs[idx] = -99
        elif inistates[idx] == 1:
            bs[idx] = 99
    inistates = np.array(inistates)
    inistates[inistates == 9] = 0
    objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
    objboolmask = objstates != 9
    objphen = objstates[objboolmask]
    
    tgtlist = []
    for drugtgt in range(1, maxdrugtgt + 1):
        contgts = it.combinations(conlist, drugtgt)
        for contgt in contgts:
            conbs = bs.copy()
            for tgt in contgt:
                conbs[node_list.index(tgt)] = sign * 99###direction of control
            attractor = _func(conbs, inistates)
            tmpfitlist = np.zeros((len(attractor),))
            for i, predstates in enumerate(attractor):
                predstates = np.array(predstates)
                predphen = predstates[objboolmask]
                tmpfit = np.sum((objphen - predphen)**2)
                tmpfitlist[i] = tmpfit
            fitness = np.mean(tmpfitlist)
            if fitness == 0:
                #check if the combination of the targets is already in the target drug list
                unique = True
                for tgtno in range(1, len(contgt)):
                    drugcombs = it.combinations(contgt, tgtno)
                    if not unique:
                        break
                    for drugcomb in drugcombs:
                        if drugcomb in tgtlist:
                            unique = False
                            break
                if unique:
                    tgtlist.append(contgt)
    tgtlists.append(tgtlist)

if mode:
    print(tgtlists[0])
else:
    def plotdata(tgtdict):
        #plot data
        plt.rcdefaults()
        fig, ax = plt.subplots()
        hbardata = []
        for node in tgtdict.keys():
            hbardata.append([tgtdict[node], node])
        hbardata.sort(reverse=True)
        if len(hbardata) > 20:
            hbardata = hbardata[:20]
        y_pos = np.arange(len(hbardata))
        ax.barh(y_pos, [num for num, _ in hbardata], align='center')
        ax.set_yticks(y_pos)
        y_labels = []
        for _, node in hbardata:
            if isinstance(node, str): 
                y_labels.append(str(node))
            else:
                y_labels.append(' & '.join([str(tgt)  for tgt in node]))
                
        ax.set_yticklabels(y_labels, fontsize=14)
        ax.invert_yaxis()
        ax.set_xlabel('Number of Occurrence', fontsize=14)
        plt.tick_params(labelsize=14)
        #ax.set_title(dataname)
        ax.set_title('Quiescence to Senescence Control: Knockout (211 Models)', fontsize=20)
        plt.tight_layout()
        plt.show()
    
    #mode of individual drug targets
    indivmode = {}
    for node in node_list:
        indivmode[node] = 0
    for tgtlist in tgtlists:
        for tgt in tgtlist:
            for node in tgt:
                indivmode[node] += 1

    #mode of tuples of drug combinations
    tuplemode = {}
    for tgtlist in tgtlists:
        for tgt in tgtlist:
            if tgt in tuplemode:
                tuplemode[tgt] += 1
            else:
                tuplemode[tgt] = 0

    #plot data
    plotdata(indivmode)
    plotdata(tuplemode)







