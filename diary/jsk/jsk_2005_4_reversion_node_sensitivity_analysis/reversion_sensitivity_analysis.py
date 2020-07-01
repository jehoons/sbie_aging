import numpy as np
import pandas as pd
import itertools as it

#algorithm parameters
sign = -1#direction of control; -1 or 1
maxdrugtgt = 1#maximum number of drug used
n = 5000

#import
netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

poisson_lambda_min = 1
poisson_lambda_max = 7
poisson_lambda_step = .25
max_disruption = 20
repeat = 20
for poisson_lambda in np.arange(poisson_lambda_min, poisson_lambda_max+poisson_lambda_step, poisson_lambda_step):
    for dataname in ['obj_attractor_control_sq.csv', 'obj_attractor_control_sp.csv']:
        conddata = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
        node_list = list(conddata.columns)
        
        iptnodes = ['lowNutrition', 'IGF1', 'DNAdamage']
        optnodes = ['S6K1', 'EIF4EBP1', 'ULK1', 'IL6', 'IL1B', 'E2F1']
        signalnodes = [node for node in node_list if not node in iptnodes and not node in optnodes]
        
        sensitivity_result = pd.DataFrame(index=signalnodes)
        for disruption in range(0, max_disruption + 1):
            for s in range(repeat):
                if disruption == 0 and s > 0:
                    continue
                print('node {} sensitivity: {}/{}, iteration: {}/{}, lambda: {}/{}'.format(dataname[-6:-4], disruption, max_disruption, s+1, repeat, float(poisson_lambda), float(poisson_lambda_max)))
                np.random.seed(s)
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
                
                def generate_fake_node(netdata_, node_name, candidate_nodes, poisson_lambda):
                    links = []
                    #choosing indegrees
                    num_indeg = 0
                    while not(num_indeg):
                        num_indeg = np.random.poisson(poisson_lambda)
                    indeg_nodes = np.random.choice(candidate_nodes, num_indeg, replace=False)
                    
                    for indeg_node in indeg_nodes:
                        interaction = np.random.choice(['activate', 'inhibit'])
                        links.append([indeg_node, interaction, node_name])
                    #choosing outdegrees
                    num_outdeg = 0
                    while not(num_outdeg):
                        num_outdeg = np.random.poisson(poisson_lambda)
                    outdeg_nodes = np.random.choice(candidate_nodes, num_outdeg, replace=False)
                    
                    for outdeg_node in outdeg_nodes:
                        interaction = np.random.choice(['activate', 'inhibit'])
                        links.append([node_name, interaction, outdeg_node])
                    return links
                
                tgtlists = []
                attractorlists = []
                statehistorylists = []
                for netnum, (ws, bs) in enumerate(netparamslist):
                    _netdata = netdata.copy()
                    _node_list = node_list.copy()
                    
                    #adding synthetic nodes
                    signalnodes_ = signalnodes.copy()
                    num_added_links = 0
                    for i in range(disruption):
                        node_name_ = 'node_{}'.format(i)
                        _links = generate_fake_node(_netdata, node_name_, signalnodes_, poisson_lambda=poisson_lambda)
                        #add links and adjust weights
                        _netdata = _netdata + _links
                        ws = np.concatenate((ws, np.random.rand(len(_links)) * 18))
                        #add node and adjust basal activity
                        _node_list.append(node_name_)
                        bs = np.append(bs, np.random.choice([1, -1]) * np.random.rand()*18)
                        #manage for loop
                        num_added_links += len(_links)
                        signalnodes_.append(node_name_)
                    
                    signs = np.ones((num_edges+num_added_links,), dtype=np.int32)
                    ind_rows = np.zeros((num_edges+num_added_links,), dtype=np.int32)
                    ind_cols = np.zeros((num_edges+num_added_links,), dtype=np.int32)
                    
                    for idx, edge in enumerate(_netdata):
                        if edge[1] == 'inhibit':
                            signs[idx] = -1
                        ind_rows[idx] = _node_list.index(edge[0])
                        ind_cols[idx] = _node_list.index(edge[2])
                    # end of for
                    
                    def _index_hist(st, hist):
                        for i, state in enumerate(hist):
                            if (state == st).all():
                                return i
                        return -1
                    
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
                    
                    inistates = conddata.loc['INPUT'].tolist()
                    for idx in iptidx:#fixing input nodes using extreme basal activity settings
                        if inistates[idx] == 0:
                            bs[idx] = -99
                        elif inistates[idx] == 1:
                            bs[idx] = 99
                    inistates = np.array(inistates)
                    inistates[inistates == 9] = 0
                    objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
                    inistates = np.concatenate((inistates, np.zeros(disruption)))
                    objstates = np.concatenate((objstates, np.array([9]*disruption)))
                    objboolmask = objstates != 9
                    objphen = objstates[objboolmask]
                    
                    tgtlist = []
                    attractorlist = []
                    statehistorylist = []
                    conbs = bs.copy()
                    attractor, _ = _func(conbs, inistates)
                    inistates = attractor[0].copy()
                    for drugtgt in range(1, maxdrugtgt + 1):
                        contgts = it.combinations(conlist, drugtgt)
                        for contgt in contgts:
                            conbs = bs.copy()
                            for tgt in contgt:
                                conbs[_node_list.index(tgt)] = sign * 99###direction of control
                            attractor, state_history = _func(conbs, inistates)
                
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
                                    attractorlist.append(attractor)
                                    tgtlist.append(contgt)
                                    statehistorylist.append(state_history)
                    tgtlists.append(tgtlist)
                    attractorlists.append(attractorlist)
                    statehistorylists.append(statehistorylist)
                
                #mode of tuples of drug combinations
                tuplemode = {}
                for tgtlist in tgtlists:
                    for tgt in tgtlist:
                        if maxdrugtgt == 1:
                            tgt = tgt[0]
                        if tgt in tuplemode:
                            tuplemode[tgt] += 1
                        else:
                            tuplemode[tgt] = 0
                
                #store sensitivity analysis data
                sensitivity_result = pd.concat([sensitivity_result, pd.Series(tuplemode, name='node_added_{}_{}'.format(disruption, s))], axis=1, sort=True).fillna(0)
        
        sensitivity_result.to_csv('result/node_sensitivity_result_{}_lambda_{}.csv'.format(dataname[-6:-4], poisson_lambda))
















