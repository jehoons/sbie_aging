import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt

#algorithm parameters
n = 5000

#import
netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

dataname = 'obj_attractor.csv'
conddata_original = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata_original.columns)

attractor_avgs = []
state_names = []
for iptcond in ['s', 'q', 'p']:
    if iptcond == 's':
        state_name = 'Senescence'
        conddata = conddata_original.iloc[[0,3]]
        conddata.index = ['INPUT', 'OUTPUT']
    elif iptcond == 'q':
        state_name = 'Quiescence'
        conddata = conddata_original.iloc[[1,4]]
        conddata.index = ['INPUT', 'OUTPUT']
    elif iptcond == 'p':
        state_name = 'Proliferation'
        conddata = conddata_original.iloc[[2,5]]
        conddata.index = ['INPUT', 'OUTPUT']
    else:
        raise ValueError
    state_names.append(state_name)
    
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
    
    ##########
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
    ##########
    
    attractorlist = []
    statehistorylist = []
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
                    return attractor, state_history
                
                else: # move on
                    state_history.append(x_t2.copy())
                    x_t1 = x_t2
        
        inistates = conddata.loc['INPUT'].tolist()
        inistates = np.array(inistates)
        inistates[inistates == 9] = 0
        objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
        objboolmask = objstates != 9
        objphen = objstates[objboolmask]
        
        attractor, state_history = _func(bs, inistates)
        tmpfitlist = np.zeros((len(attractor),))
        for i, predstates in enumerate(attractor):
            predstates = np.array(predstates)
            predphen = predstates[objboolmask]
            tmpfit = np.sum((objphen - predphen)**2)
            tmpfitlist[i] = tmpfit
        fitness = np.mean(tmpfitlist)
        
        assert fitness == 0
        
        attractorlist.append(attractor)
        statehistorylist.append(state_history)
    
    
    #tracing input signal propagation
    state_len_hist_plot = False
    fitness_by_time_step_plot = False
    
    #plot signal propagation time steps
    statehislens = np.array([len(statehis) for statehis in statehistorylist])
    if state_len_hist_plot:
        plt.figure()
        plt.title('{}: Signal Propagation Time Steps'.format(state_name))
        plt.xlabel('time steps')
        plt.ylabel('models')
        plt.hist(statehislens)
    
    #plot fitness against time step of the ensemble model
    fithistorylist = []
    maxshl = max(statehislens)
    for state_history in statehistorylist:
        fit_history = []
        for j, state in enumerate(state_history):
            statephen = np.array(state)[objboolmask]
            interfit = np.sum((objphen - statephen)**2).astype(np.float64)
            fit_history.append(interfit)
        j += 1
        while j < maxshl:
            fit_history.append(interfit)
            j += 1
        fithistorylist.append(fit_history)
    meanfitbyts = np.mean(np.array(fithistorylist), axis=0)
    fitstdbyts = np.std(np.array(fithistorylist), axis=0)
    fitupper = meanfitbyts + fitstdbyts
    fitlower = np.maximum(meanfitbyts - fitstdbyts, 0)
    
    if fitness_by_time_step_plot:
        plt.figure()
        plt.title('{}: Fitness by Time Steps'.format(state_name))
        plt.xlabel('time steps')
        plt.ylabel('fitness')
        timesteps = np.arange(1,1+len(meanfitbyts))
        plt.plot(timesteps, meanfitbyts, color='b', lw=2, alpha=.8)
        plt.fill_between(timesteps, fitlower, fitupper, color='grey', alpha=.2)
    
    #attractor average of all networks in the ensemble
    _attractor = []
    for attractors in attractorlist:
        _attractor.append(np.mean(attractors, axis=0))
    attractor_avg = np.mean(_attractor, axis=0)
    #print(attractor_avg)
    attractor_avgs.append(list(attractor_avg))
attractor_avgs = np.array(attractor_avgs)

pd.DataFrame(attractor_avgs, index=state_names, columns=node_list).to_csv('attractor_avgs.csv')















































