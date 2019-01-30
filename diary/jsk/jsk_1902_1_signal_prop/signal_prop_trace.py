import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle


plot = False

#import data
with open('target_list.bin', 'rb') as f:
    targetlist = pickle.load(f)
with open('attractor_list.bin', 'rb') as f:
    attractorlist = pickle.load(f)
with open('state_history_list.bin', 'rb') as f:
    statehistorylist = pickle.load(f)


#mark PDK1 controlled models
pdk1attractors = []
pdk1statehistory = []
for idx, targets in enumerate(targetlist):
    if ('PDK1',) in targets:
        pdk1attractors.append(attractorlist[idx][targets.index(('PDK1',))])
        pdk1statehistory.append(statehistorylist[idx][targets.index(('PDK1',))])

statehislens = np.array([len(statehis) for statehis in pdk1statehistory])
if plot:
    plt.hist(statehislens)

#tentatively decided to only include statehistories with length 4-8
mask = np.logical_and(statehislens >= 4, statehislens <= 8)
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
for stlen in range(4, 9):
    state_hist[stlen] = state_histories[state_hist_len == stlen]
    indiv_mst[stlen] = mean_state_transition(state_hist[stlen])

#altogether individually
#fill empty ends with attractor states
state_hist_cum = []
for stlen in range(4, 9):
    for _states in state_hist[stlen]:
        state_hist_cum.append(_states.copy() + [_states[-1]]*(8-stlen))
cum_mst = mean_state_transition(state_hist_cum)

#normalized by steps
binsize = np.lcm.reduce(np.arange(3,8))
bins = np.zeros((binsize+1, attlen))
state_hist_norm = []
for _states in state_histories:
    interpol = bins.copy()
    minibinsize = int(binsize/(len(_states)-1))
    for i in range(len(_states) - 1):
        for j in range(minibinsize):
            interpol[i*minibinsize+j] = (_states[i+1]-_states[i])/minibinsize*j + _states[i]
    #handle obob
    interpol[-1] = _states[-1]
    
    _norm_states = []
    for k in range(8):
        _norm_states.append(interpol[int(binsize/7*k)])
    state_hist_norm.append(_norm_states)
norm_mst = mean_state_transition(state_hist_norm)

state_transition_sigs = {'individually': indiv_mst, 'cumulative': cum_mst, 'normalized': norm_mst}

#save
with open('selected_state_histories.bin', 'wb') as f:
    pickle.dump(state_histories, f)
with open('state_transition_signatures.bin', 'wb') as f:
    pickle.dump(state_transition_sigs, f)
















