import numpy as np
import pandas as pd
import pickle


"""
remove positive feedback for logically removing PTEN's signal flow
"""

'''
from pyeda import *

bool_dict = {}
for boolean_output in ["{0:b}".format(i).zfill(4) for i in range(2**4)]:
    x = exprvars('x', 2)
    f = truthtable(x, boolean_output)
    logic = espresso_tts(f)
    #logic = truthtable2expr(f)
    bool_dict[boolean_output] = [str(logic)]
bool_df = pd.DataFrame.from_dict(bool_dict).T
bool_df.to_csv('bool.csv', sep='\t', header=False)
'''

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

#for pdk1_controlled_network_only in [False, True]:
for positive_feedback_network_only in [False, True]:
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
    
#    if pdk1_controlled_network_only:
#        #import data
#        with open('target_list.bin', 'rb') as f:
#            targetlist = pickle.load(f)
#        
#        #mark PDK1 controlled models
#        pdk1netparamslist = []
#        for idx, targets in enumerate(targetlist):
#            if ('PDK1',) in targets:
#                pdk1netparamslist.append(netparamslist[idx])
#        netparamslist = pdk1netparamslist
    
    if positive_feedback_network_only:
        #import data
        with open('positive_feedback_model.bin', 'rb') as f:
            pos_feed_net_idxs = pickle.load(f)
        pos_feed_netparamslist = []
        for idx in pos_feed_net_idxs:
            pos_feed_netparamslist.append(netparamslist[idx])
        netparamslist = pos_feed_netparamslist
    
    #run Boolean network simulation
    def _index_hist(st, hist):
        for i, state in enumerate(hist):
            if (state == st).all():
                return i
        return -1
    
    #weighted-sum to logic transformation for pdk1
    pik3ca_idx = node_list.index('PIK3CA')
    pdk1_idx = node_list.index('PDK1')
    pik3ca_pdk1_widx = netdata.index(['PIK3CA', 'activate', 'PDK1'])
    pten_pdk1_widx = netdata.index(['PTEN', 'inhibit', 'PDK1'])
    bool_strs = []
    for ws, bs in netparamslist:
        ws = np.round(ws)
        bool_str = ''
        for pik3ca, pten in [(0,0), (0,1), (1,0), (1,1)]:
            pdk1 = ws[pik3ca_pdk1_widx] * pik3ca - ws[pten_pdk1_widx] * pten + bs[pdk1_idx]
            if pdk1 > 0:
                bool_str = bool_str + '1'
            else:
                bool_str = bool_str + '0'
        bool_strs.append(bool_str)
    bool_strs_df = pd.DataFrame.from_dict({'pik3ca, pten: (00, 01, 10, 11)': bool_strs})
    
    
    pdk1_logic_dict = {'0000': lambda pdk1: 0, '0010': lambda pdk1: pdk1, '0011': lambda pdk1: pdk1,
                       '1010': lambda pdk1: -1, '1011': lambda pdk1: pdk1, '1111': lambda pdk1: 1}

    for pik3ca_inhibit, pos_feed_cut in zip([False, False, True, True], [False, True, False, True]):
        attractorlist = []
        statehistorylist = []
        for (ws, bs), bool_str in zip(netparamslist, bool_strs):
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
            def _func1(b, inistates):
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
            
            attractor, statehistory = _func1(bs, inistates)
            attractorlist.append(attractor)
            statehistorylist.append(statehistory)
        
        #analyze resulting attractor information
        attractor_df = pd.DataFrame(index=np.arange(len(attractorlist)), columns=node_list)
        for _i, _attractor in enumerate(attractorlist):
            attractor_df.loc[_i] = np.mean(_attractor, axis=0)
        attractor_avg = attractor_df.mean()
        positive_feedback_bf = attractor_avg.loc[positive_feedback_nodes]
        
        ipt_condition_new_list = []
        attractor_df_new_list = []
        attractor_avg_new_list = []
        positive_feedback_af_list = []
        pos_feed_link_idx1 = netdata.index(['AKT1', 'inhibit', 'PTEN'])
        pos_feed_link_idx2 = netdata.index(['IKBKB', 'inhibit', 'PTEN'])
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
                def _func2(w, b, inistates):
                    weight_matrix = np.zeros((len(_node_list), len(_node_list)))
                    weight_matrix[ind_rows, ind_cols] = signs*np.round(w)
                    
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
                
                        # Force positive feedback cut; PTEN -| PDK1
                        if pos_feed_cut:
                            x_t2[pdk1_idx] = pdk1_logic_dict[bool_str](x_t1[pik3ca_idx])
                
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
                
                #force pik3ca
                if pik3ca_inhibit:
                    inistates[pik3ca_idx] = 0
                    bs[pik3ca_idx] = -99
                
                attractor, statehistory = _func2(ws, bs, inistates)
                attractorlist_new.append(attractor)
                statehistorylist_new.append(statehistory)
            
            #analyze resulting attractor information
            attractor_df_new = pd.DataFrame(index=np.arange(len(attractorlist_new)), columns=node_list)
            for _i, _attractor in enumerate(attractorlist_new):
                attractor_df_new.loc[_i] = np.mean(_attractor, axis=0)
            attractor_avg_new = attractor_df_new.mean()
            positive_feedback_af = attractor_avg_new.loc[positive_feedback_nodes]
            
            ipt_condition_new_list.append(pd.Series([int(_cond) for _cond in _ipt_cond], index=np.array(node_list)[iptidx]))
            attractor_df_new_list.append(attractor_df_new)
            attractor_avg_new_list.append(attractor_avg_new)
            positive_feedback_af_list.append(positive_feedback_af)
        
        #organize results
        ipt_condition_new_df = pd.concat(ipt_condition_new_list, axis=1)
        attractor_avg_new_df = pd.concat(attractor_avg_new_list, axis=1)
        positive_feedback_af_df = attractor_avg_new_df.loc[positive_feedback_nodes]
        output_nodes_af_df = attractor_avg_new_df.loc[np.array(node_list)[optidx]]
        
        all_df = pd.concat([ipt_condition_new_df, positive_feedback_af_df, output_nodes_af_df])
        
        if pos_feed_cut:
            pos_feed_cut_str = 'cut'
        else:
            pos_feed_cut_str = 'intact'
        if pik3ca_inhibit:
            pik3ca_inhibit_str= '_pik3ca_inhibited'
        else:
            pik3ca_inhibit_str = ''
#        if pdk1_controlled_network_only:
#            netparams_str = '_pdk1_network_only'
#        else:
#            netparams_str = ''
        if positive_feedback_network_only:
            netparams_str = '_pos_feed_network_only'
        else:
            netparams_str = ''
        
        #save organized and total result
        all_df.to_csv('results_bool/simulation_result_pos_feed_{}{}{}.csv'.format(pos_feed_cut_str, pik3ca_inhibit_str, netparams_str))
        attractor_avg_new_df.to_csv('results_bool/average_attractor_pos_feed_{}{}{}.csv'.format(pos_feed_cut_str, pik3ca_inhibit_str, netparams_str))
        
        iptnodes = np.array(node_list)[iptidx]
        for _ipt_cond, _attractor_df_new in zip(["{0:b}".format(i).zfill(3) for i in range(2**3)], attractor_df_new_list):
            _attractor_df_new.to_csv('results_bool/attractor_res_raw_{}_{}_{}_{}_{}_{}_{}{}{}.csv'.format(iptnodes[0], _ipt_cond[0], iptnodes[1], _ipt_cond[1], iptnodes[2], _ipt_cond[2], pos_feed_cut_str, pik3ca_inhibit_str, netparams_str))
    bool_strs_df.to_csv('results_bool/logic_conversion_table{}.csv'.format(netparams_str))











