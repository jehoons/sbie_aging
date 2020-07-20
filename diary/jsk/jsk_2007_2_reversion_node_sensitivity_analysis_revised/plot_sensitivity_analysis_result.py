import numpy as np
import pandas as pd
from scipy.stats import weightedtau
import matplotlib.pyplot as plt
import seaborn as sns


#parameters
max_int = 20
rank_plot_threshold = 8
plot_boxplot = False
plot_lineplot = True
lambdas_2_plot = [2, 4, 6]
plot_rank_movement = False
plot_3d = True
plot_save = True

#weigher function for weighted tau rank correlation calculation
weigher = lambda x: 1/(x+1)**2

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
num_node = len(node_list)


poisson_lambda_min = 1
poisson_lambda_max = 7
poisson_lambda_step = .25
data = {}
dfs = []
for poisson_lambda in np.arange(poisson_lambda_min, poisson_lambda_max+poisson_lambda_step, poisson_lambda_step):
    data[poisson_lambda] = {}
    for sq_or_sp in ['sq', 'sp']:
        data[poisson_lambda][sq_or_sp] = pd.read_csv('result/node_sensitivity_result_{}_lambda_{}.csv'.format(sq_or_sp, poisson_lambda), index_col=0)
    data[poisson_lambda]['score'] = data[poisson_lambda]['sq'] - data[poisson_lambda]['sp']
    _data = data[poisson_lambda]['score']
    
    #plot boxplots
    bpdict = {}
    for simname in _data.columns:
        _type, _, _int, _iter = simname.split('_')
        if int(_int) == max_int + 1:
            break
        if int(_int) == 0:
            baseline = _data[simname]
            continue
        if not int(_int) in bpdict.keys():
            bpdict[int(_int)] = []
        _corr = weightedtau(baseline, _data[simname], weigher=weigher)
        bpdict[int(_int)].append(_corr[0])
    
    data[poisson_lambda]['correlation'] = bpdict
    bplabels, bpdata = list(bpdict.keys()), bpdict.values()
    
    width = 1 / num_node * .8
    max_range = len(bplabels) / num_node + width
    
    if plot_boxplot:
        plt.figure()
        plt.boxplot(bpdata, positions=np.array(bplabels)/num_node, widths=width, manage_xticks=False)#widths and ticks
        plt.xticks(np.arange(0, 10*max_range+1)/10)
        plt.xlim([0, max_range])
        plt.xlabel('number of nodes added / total number of nodes')
        plt.ylabel('weighted kendall tau rank correlation')
        if plot_save:
            plt.savefig('node_boxplot.jpg')
    
    if plot_lineplot:
        sns.set(style='whitegrid')
        
        _df = pd.DataFrame(index=['fraction of added nodes [%]', 'rank correlation of identified targets'])
        _ = 0
        for _nodes_added, _kendall_taus in bpdict.items():
            _node_fraction = _nodes_added / num_node * 100
            for _kendall_tau in _kendall_taus:
                _df[_] = [_node_fraction, _kendall_tau]
                _ += 1
        df = _df.T
        
        if poisson_lambda in lambdas_2_plot:
            df['degree parameter of\nadded nodes (lambda)'] = [poisson_lambda]*df.shape[0]
            dfs.append(df)
    
    #rank node results
    _ranked = pd.DataFrame(index=np.arange(1, _data.shape[0]+1))
    _node_list = list(_data.index)
    for _sim in _data.columns:
        _ranked[_sim] = sorted(_node_list, key=lambda _node: _data[_sim][_node], reverse=True)
    data[poisson_lambda]['ranked'] = _ranked
    
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
    _rank_movement = pd.DataFrame(index=_data.index, columns=np.arange(max_int+1))
    for _int in _rank_movement.columns:
        if _int == 0:
            _data[_data.columns[0]]
        for _node in _rank_movement.index:
            _rank_movement[_int][_node] = np.mean(_tmp[_int][_node])
    _rank_movement = _rank_movement.sort_values(by=[0])
    data[poisson_lambda]['rank_by_int'] = _rank_movement
    
    if plot_rank_movement:
        plt.figure()
        plt.ylim([0, rank_plot_threshold+1])
        plt.yticks(np.arange(1,rank_plot_threshold+1))
        plt.gca().invert_yaxis()
        for _node, _ranks in _rank_movement.iterrows():
            if all(_ranks > rank_plot_threshold):
                continue
            plt.plot(_ranks, label=_node)
        plt.xticks(np.arange(0, max_int+1, 2))
        plt.legend(loc='lower right')
        plt.xlabel('number of interventions')
        plt.ylabel('node average rank')
        if plot_save:
            plt.savefig('node_rank_movement.jpg', dpi=300)

if plot_3d:
    from mpl_toolkits import mplot3d
    
    _lambdas = np.linspace(poisson_lambda_min, poisson_lambda_max, (poisson_lambda_max - poisson_lambda_min) / poisson_lambda_step + 1)
    _node_interventions = np.linspace(1, max_int, max_int)
    _node_fractions = _node_interventions / num_node * 100
    
    lambdas, node_fractions = np.meshgrid(_lambdas, _node_fractions)
    kendalltaus = np.zeros((len(_node_interventions), len(_lambdas)))
    for _x, _node_intervention in enumerate(_node_interventions):
        for _y, _lambda in enumerate(_lambdas):
            kendalltaus[_x][_y] = np.median(data[_lambda]['correlation'][_node_intervention])
    
    fig = plt.figure()
    ax=plt.axes(projection='3d')
    ax.plot_surface(lambdas, node_fractions, kendalltaus, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel('degree parameter of added nodes (lambda)')
    ax.set_ylabel('fraction of added nodes [%]')
    ax.set_zlabel('rank correlation of identified targets')
    
    if plot_save:
        plt.savefig('3d_projection.jpg', dpi=300)

if plot_lineplot:
    fig = plt.figure()
    g = sns.lineplot(x='fraction of added nodes [%]', y='rank correlation of identified targets', hue='degree parameter of\nadded nodes (lambda)', palette=sns.color_palette('deep', len(lambdas_2_plot)), ci='sd', data=pd.concat(dfs, ignore_index=True), legend='full')
    
    #bolden lambda 2.0
    l = g.lines[0]
    plt.setp(l, linewidth=2.5)
    s = g.collections[0]
    plt.setp(s, alpha=.5)
    
    if len(lambdas_2_plot) > 0:
        for i in range(1, len(lambdas_2_plot)):
            g.lines[i].set_linestyle('--')
            plt.setp(g.collections[i], alpha=.25)
    
    for leg_line in g.legend().get_lines()[2:]:
        leg_line.set_linestyle('--')
    
    if plot_save:
        fig.savefig('lineplot.jpg', dpi=300)

































