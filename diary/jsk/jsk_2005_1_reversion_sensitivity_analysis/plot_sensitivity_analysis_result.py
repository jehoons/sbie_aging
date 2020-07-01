import numpy as np
import pandas as pd
from scipy.stats import weightedtau
import matplotlib.pyplot as plt


#parameters
valid_num = 8
data_dir = 'data_lambda_{}/'.format(2.125)#1.878 or 2.125
max_int = {'edge': 20, 'node': 12}
rank_plot_threshold = 8
plt_save = True

weigher = lambda x: 1/(x+1)**2
#weigher = lambda x: np.e**(-x)

#import network data
netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))
_node_list = []
for link in netdata:
    if not link[0] in _node_list:
        _node_list.append(link[0])
    if not link[2] in _node_list:
        _node_list.append(link[2])
node_list = sorted(_node_list)

num_edge = len(netdata)
num_node = len(node_list)
num_tot = {'edge': num_edge, 'node': num_node}

data = {}
for edge_or_node in ['edge', 'node']:
    data[edge_or_node] = {}
    for sq_or_sp in ['sq', 'sp']:
        data[edge_or_node][sq_or_sp] = pd.read_csv(data_dir + 'sensitivity_result_{}_{}.csv'.format(edge_or_node, sq_or_sp), index_col=0)
    data[edge_or_node]['score'] = data[edge_or_node]['sq'] - data[edge_or_node]['sp']
    _data = data[edge_or_node]['score']
    
    #plot boxplots
    bpdict = {}
    for simname in _data.columns:
        _type, _, _int, _iter = simname.split('_')
        if int(_int) == max_int[edge_or_node] + 1:
            break
        if int(_int) == 0:
            baseline = _data[simname]
            continue
        if not int(_int) in bpdict.keys():
            bpdict[int(_int)] = []
        _corr = weightedtau(baseline, _data[simname], weigher=weigher)
        bpdict[int(_int)].append(_corr[0])
    data[edge_or_node]['correlations'] = bpdict
    
    bplabels, bpdata = list(bpdict.keys()), bpdict.values()
    
    width = 1 / num_tot[edge_or_node] * .8
    max_range = len(bplabels) / num_tot[edge_or_node] + width
    
    plt.figure()
    plt.boxplot(bpdata, positions=np.array(bplabels)/num_tot[edge_or_node], widths=width, manage_xticks=False)#widths nad ticks
    
    plt.xticks(np.arange(0, 10*max_range+1)/10)
    plt.xlim([0, max_range])
    plt.xlabel('number of {}s added / total number of {}s'.format(edge_or_node, edge_or_node))
    plt.ylabel('weighted kendall tau rank correlation')
    if plt_save:
        plt.savefig('boxplot_{}.jpg'.format(edge_or_node))
    
    #rank node results
    _ranked = pd.DataFrame(index=np.arange(1, _data.shape[0]+1))
    _node_list = list(_data.index)
    for _sim in _data.columns:
        _ranked[_sim] = sorted(_node_list, key=lambda _node: _data[_sim][_node], reverse=True)
    data[edge_or_node]['ranked'] = _ranked
    
    _tmp = {}
    for simname in _data.columns:
        _type, _, _int, _iter = simname.split('_')
        if not int(_int) in _tmp.keys():
            _tmp[int(_int)] = {}
            for _node in _data.index:
                _tmp[int(_int)][_node] = []
        for _rank, _node in _ranked[simname].items():
            _tmp[int(_int)][_node].append(_rank)
    
    #plot rank movement by intervention
    _rank_movement = pd.DataFrame(index=_data.index, columns=np.arange(max_int[edge_or_node]+1))
    for _int in _rank_movement.columns:
        if _int == 0:
            _data[_data.columns[0]]
        for _node in _rank_movement.index:
            _rank_movement[_int][_node] = np.mean(_tmp[_int][_node])
    _rank_movement = _rank_movement.sort_values(by=[0])
    data[edge_or_node]['rank_by_int'] = _rank_movement
    
    plt.figure()
    plt.ylim([0, rank_plot_threshold+1])
    plt.yticks(np.arange(1,rank_plot_threshold+1))
    plt.gca().invert_yaxis()
    for _node, _ranks in _rank_movement.iterrows():
        if all(_ranks > rank_plot_threshold):
            continue
        plt.plot(_ranks, label=_node)
    plt.xticks(np.arange(0, max_int[edge_or_node]+1, 2))
    plt.legend(loc='lower right')
    plt.xlabel('number of interventions')
    plt.ylabel('node average rank')
    if plt_save:
        plt.savefig('rank_movement_{}.jpg'.format(edge_or_node))
    
    




























