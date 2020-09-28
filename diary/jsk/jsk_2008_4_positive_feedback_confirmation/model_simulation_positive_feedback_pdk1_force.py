import numpy as np
import pandas as pd
import pickle


#algorithm parameters
n = 5000
#positive_feedback = {'PDK1': 1, 'AKT': 1, 'IKBKB': 1, 'PTEN': 0}
positive_feedback_nodes = ['PDK1', 'AKT1', 'IKBKB', 'PTEN']
#pdk1_inhibit = True
#pdk1_controlled_network_only = False

for pdk1_inhibit, pdk1_controlled_network_only in zip([False, True, False, True], [False, False, True, True]):
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
    
    if pdk1_controlled_network_only:
        #import data
        with open('target_list.bin', 'rb') as f:
            targetlist = pickle.load(f)
        
        #mark PDK1 controlled models
        pdk1netparamslist = []
        for idx, targets in enumerate(targetlist):
            if ('PDK1',) in targets:
                pdk1netparamslist.append(netparamslist[idx])
        netparamslist = pdk1netparamslist
    
    
    #run Boolean network simulation
    def _index_hist(st, hist):
        for i, state in enumerate(hist):
            if (state == st).all():
                return i
        return -1
    
    
    pdk1_idx = node_list.index('PDK1')
    for pdk1_force in [(False, False), (False, True), (True, False), (True, True)]:
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
            
            #force pdk1
            if pdk1_force[0]:
                if pdk1_inhibit:
                    inistates[pdk1_idx] = 0
                    bs[pdk1_idx] = -99
                else:
                    inistates[pdk1_idx] = 1
                    bs[pdk1_idx] = 99
            
            attractor, statehistory = _func(bs, inistates)
            attractorlist.append(attractor)
            statehistorylist.append(statehistory)
        
        
        #analyze resulting attractor information
        attractor_df = pd.DataFrame(index=np.arange(len(attractorlist)), columns=node_list)
        for _i, _attractor in enumerate(attractorlist):
            attractor_df.loc[_i] = np.mean(_attractor, axis=0)
        attractor_avg = attractor_df.mean()
        positive_feedback_bf = attractor_avg.loc[positive_feedback_nodes]
        
        ipt_condition_new_list = []
        attractor_avg_new_list = []
        positive_feedback_af_list = []
        #rerun senescence attractors for all input conditions
        for _ipt_cond in ["{0:b}".format(i).zfill(3) for i in range(2**3)]:
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
                for _idx, _cond in zip(iptidx, _ipt_cond):#fixing input nodes using extreme basal activity settings
                    inistates[_idx] = int(_cond)
                    if int(_cond):
                        bs[_idx] = 99
                    else:
                        bs[_idx] = -99
                
                #force pdk1
                if pdk1_force[1]:
                    if pdk1_inhibit:
                        inistates[pdk1_idx] = 0
                        bs[pdk1_idx] = -99
                    else:
                        inistates[pdk1_idx] = 1
                        bs[pdk1_idx] = 99
                
                attractor, statehistory = _func(bs, inistates)
                attractorlist_new.append(attractor)
                statehistorylist_new.append(statehistory)
            
            #analyze resulting attractor information
            attractor_df_new = pd.DataFrame(index=np.arange(len(attractorlist_new)), columns=node_list)
            for _i, _attractor in enumerate(attractorlist_new):
                attractor_df_new.loc[_i] = np.mean(_attractor, axis=0)
            attractor_avg_new = attractor_df_new.mean()
            positive_feedback_af = attractor_avg_new.loc[positive_feedback_nodes]
            
            ipt_condition_new_list.append(pd.Series([int(_cond) for _cond in _ipt_cond], index=np.array(node_list)[iptidx]))
            attractor_avg_new_list.append(attractor_avg_new)
            positive_feedback_af_list.append(positive_feedback_af)
        
        #organize results
        ipt_condition_new_df = pd.concat(ipt_condition_new_list, axis=1)
        attractor_avg_new_df = pd.concat(attractor_avg_new_list, axis=1)
        positive_feedback_af_df = attractor_avg_new_df.loc[positive_feedback_nodes]
        output_nodes_af_df = attractor_avg_new_df.loc[np.array(node_list)[optidx]]
        
        all_df = pd.concat([ipt_condition_new_df, positive_feedback_af_df, output_nodes_af_df])
        if pdk1_inhibit:
            pdk1_force_str = 'inhibit'
        else:
            pdk1_force_str = 'activate'
        if pdk1_controlled_network_only:
            netparams_str = '_pdk1_network_only'
        else:
            netparams_str = ''
        
        #save organized and total result
        all_df.to_csv('results/simulation_result_pdk1_{}_bf_{}_af_{}{}.csv'.format(pdk1_force_str, int(pdk1_force[0]), int(pdk1_force[1]), netparams_str))
        attractor_avg_new_df.to_csv('results/average_attractor_pdk1_{}_bf_{}_af_{}{}.csv'.format(pdk1_force_str, int(pdk1_force[0]), int(pdk1_force[1]), netparams_str))


































