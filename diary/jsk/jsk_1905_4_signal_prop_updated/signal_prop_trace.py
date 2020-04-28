import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle


state_len_hist_plot = True
fitness_by_time_step_plot = True
n = 5000

#import data
with open('target_list.bin', 'rb') as f:
    targetlist = pickle.load(f)
with open('attractor_list.bin', 'rb') as f:
    attractorlist = pickle.load(f)
with open('state_history_list.bin', 'rb') as f:
    statehistorylists = pickle.load(f)

netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

dataname = 'obj_attractor_control_sq.csv'
conddata = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)
objstates = np.array(conddata.loc['OUTPUT'].tolist())#objective states
objboolmask = objstates != 9
objphen = objstates[objboolmask]

#analysis on general (mean) dynamics of the network
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


#mark PDK1 controlled models
pdk1attractors = []
pdk1statehistory = []
pdk1netparamslist = []
for idx, targets in enumerate(targetlist):
    if ('PDK1',) in targets:
        pdk1attractors.append(attractorlist[idx][targets.index(('PDK1',))])
        pdk1statehistory.append(statehistorylists[idx][targets.index(('PDK1',))])
        pdk1netparamslist.append(netparamslist[idx])

statehislens = np.array([len(statehis) for statehis in pdk1statehistory])
minstep = min(statehislens)
maxstep = max(statehislens)
if state_len_hist_plot:
    plt.figure()
    plt.title('PDK1 Inhibition: Signal Propagation Time Steps')
    plt.xlabel('time steps')
    plt.ylabel('models')
    plt.hist(statehislens)
    plt.xlim([minstep-.5, maxstep+.5])

#plot fitness against time step of the ensemble model
fithistorylist = []
maxshl = max(statehislens)
for state_history in pdk1statehistory:
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
    plt.title('PDK1 Inhibition: Fitness by Time Steps')
    plt.xlabel('time steps')
    plt.ylabel('fitness')
    timesteps = np.arange(1,1+len(meanfitbyts))
    plt.plot(timesteps, meanfitbyts, color='b', lw=2, alpha=.8)
    plt.fill_between(timesteps, fitlower, fitupper, color='grey', alpha=.2)
'''
#tentatively decided to only include statehistories with length 6-10
mask = np.logical_and(statehislens >= minstep, statehislens <= maxstep)
state_histories = np.array(pdk1statehistory)[mask]

def mean_state_transition(state_library):
    #need state_transition length to be equal
    library_len = len(state_library[0])
    attlen = len(state_library[0][0])
    mst = np.zeros((library_len, attlen))
    for _state_hist in state_library:
        for idx, state in enumerate(_state_hist):
            mst[idx] += state
    mst /= len(state_library)
    return mst

attlen = len(state_histories[0][0])

#individually
state_hist_len = np.array([len(statehis) for statehis in state_histories])
state_hist = {}
indiv_mst = {}
for stlen in range(minstep, maxstep + 1):
    if not stlen in state_hist_len:
        continue
    state_hist[stlen] = state_histories[state_hist_len == stlen]
    indiv_mst[stlen] = mean_state_transition(state_hist[stlen])

#altogether individually
#fill empty ends with attractor states
state_hist_cum = []
for stlen in state_hist:
    for _states in state_hist[stlen]:
        state_hist_cum.append(_states.copy() + [_states[-1]]*(maxstep-stlen))
cum_mst = mean_state_transition(state_hist_cum)

#normalized by steps
binsize = np.lcm.reduce(np.arange(minstep, maxstep+1))
bins = np.zeros((binsize+1, attlen))
state_hist_norm = []
for _states in state_histories:
    interpol = bins.copy()
    minibinsize = int(binsize/(len(_states)-1))
    for i in range(len(_states) - 1):
        for j in range(minibinsize):
            interpol[i*minibinsize+j] = (_states[i+1]-_states[i]) / minibinsize * j + _states[i]
    #handle obob
    interpol[-1] = _states[-1]
    
    _norm_states = []
    for k in range(maxstep):
        _norm_states.append(interpol[int(binsize/(maxstep-1)*k)])
    state_hist_norm.append(_norm_states)
norm_mst = mean_state_transition(state_hist_norm)

state_transition_sigs = {'individually': indiv_mst, 'cumulative': cum_mst, 'normalized': norm_mst}

#save
with open('selected_state_histories.bin', 'wb') as f:
    pickle.dump(state_histories, f)
with open('state_transition_signatures.bin', 'wb') as f:
    pickle.dump(state_transition_sigs, f)
'''

#attractor average of all networks in the ensemble
_attractor = []
for attractors in pdk1attractors:
    _attractor.append(np.mean(attractors, axis=0))
pdk1attractor_avg = np.mean(_attractor, axis=0)

#average signal magnitude at attractor state
avgws = np.mean(np.array(pdk1netparamslist)[:,0])
signal_strength = []
for link, w in zip(netdata, avgws):
    if link[1] == 'inhibit':
        w = -w
    state = 2 * pdk1attractor_avg[node_list.index(link[0])] - 1
    signal_strength.append(state * w)
ss_comprehensive = np.hstack((np.array(netdata),np.array(signal_strength).reshape(len(netdata),1)))

pd.DataFrame(pdk1attractor_avg, index=node_list, columns=['PDK1 Inhibition']).to_csv('pdk1_inihibition_attractor_avgs.csv')
pd.DataFrame(ss_comprehensive, columns=['input', 'interaction', 'output', 'PDK1 Inhibition']).to_csv('pdk1_inihibition_signal_strength.csv')









