import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import pickle
import scipy.cluster.hierarchy as sch
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import TSNE
from matplotlib.colors import ListedColormap

'''
network_total = 5000

#import
netfilename = "PKN_24.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

with open('target_list.bin', 'rb') as f:
    targetlist = pickle.load(f)
with open('attractor_list.bin', 'rb') as f:
    attractorlist_drugged = pickle.load(f)
with open('network_parameter_list.bin', 'rb') as f:
    netparamslist = pickle.load(f)


dataname = 'obj_attractor.csv'
conddata = pd.read_csv(dataname, index_col = 0)#pd.DataFrame with input and output condition information
node_list = list(conddata.columns)

iptidx = []#indicies of the input conditions
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)

num_nodes = len(node_list)
num_edges = len(netdata)

#calculating maximum indegree to set the bounds below
indeglist = {}
for node in node_list:
    indeglist[node] = 0
for link in netdata:
    indeglist[link[2]] += 1
maxindeg = indeglist[max(indeglist, key = indeglist.get)]
wr = 2 * (maxindeg + 2)


###############################################################################
###the Boolean network simulator
signs = np.ones((num_edges,), dtype=np.int32)

ind_rows = np.zeros((num_edges,), dtype=np.int32)
ind_cols = np.zeros((num_edges,), dtype=np.int32)

for idx, edge in enumerate(netdata):
    if edge[1] == 'inhibit':
        signs[idx] = -1
    ind_rows[idx] = node_list.index(edge[0])
    ind_cols[idx] = node_list.index(edge[2])
# end of for
weight_matrix = np.zeros((num_nodes, num_nodes))


#the Boolean network model
def _index_hist(st, hist):
    for i, state in enumerate(hist):
        if (state == st).all():
            return i
    return -1

def _func(w, b, inistates):
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


attractorlist = []
incornodeslist = []
for iws, ibs in netparamslist:
    attractors = []
    incornodes = []
    for condno in range(int(len(conddata)/2)):
        condno += 1
        inistates = conddata.loc['INPUT: ' + str(condno)].tolist()
        for idx in iptidx:#fixing input nodes using extreme basal activity settings
            if inistates[idx] == 0:
                ibs[idx] = -wr
            elif inistates[idx] == 1:
                ibs[idx] = wr
        inistates = np.array(inistates)
        inistates[inistates == 9] = 0
        attractor = _func(iws, ibs, inistates)#calculated states form the Boolean model
        attractors.append(attractor)
        objstates = conddata.loc['OUTPUT: ' + str(condno)].tolist()#objective states
        objidx = [idx for idx, state in enumerate(objstates) if state != 9]
        incornode = []
        for attstate in attractor:
            for idx in objidx:
                if objstates[idx] != attstate[idx]:
                    incornode.append(node_list[idx])
            incornode.append('\t')
        incornodes.append(incornode)
    attractorlist.append(attractors)#attractors for best fit parameters
    incornodeslist.append(incornodes)#incorrect nodes for best fit parameters; none if perfect

attlist = []
for attractors in attractorlist:
    att = np.zeros((3, num_nodes))
    for i, attractor in enumerate(attractors):
        att[i] = np.mean(attractor, axis=0)
    att = att.flatten()
    attlist.append(att)
attlist = np.array(attlist)#[:network_total]


###############################################################################
###clustering analysis-attractors

#clustering conditions
save = True
clust_method = 'mandist'#'corr' or 'mandist'
plot_tsne = True
tsne_dim = 3#2 or 3 allowed

#attractor data
data = pd.DataFrame(attlist).T

#generate correlation matrix
corrmat = pd.DataFrame(0, index = data.columns, columns = data.columns).astype('float64')
for gene1 in data.columns:
    for gene2 in data.columns:
        if clust_method == 'corr':
            corrmat[gene1][gene2] = data[gene1].corr(data[gene2])
        elif clust_method == 'mandist':
            corrmat[gene1][gene2] = sum(data[gene1] - data[gene2])
        else:
            raise ValueError

#compute dendrogram
D=corrmat.values
D=np.nan_to_num(D)
Y = sch.linkage(D, method='centroid')
'''
for condition in range(6):#'Senescence': 0, 'Quiescence': 1, 'Proliferation': 2
    #the plots are run again for pdk1 marked plots
    if condition >= 3:
        mark_pdk1 = True
        condition -= 3
    else:
        mark_pdk1 = False
    
    #plot dendrogram
    fig = pylab.figure()
    axdendro = fig.add_axes([0.09,0.1,0.1,0.75])
    if mark_pdk1:
        Z = sch.dendrogram(Y, orientation='left', color_threshold=None, link_color_func=lambda x: link_colors[x])
    else:
        Z = sch.dendrogram(Y, orientation='left', color_threshold=.9*max(Y[:,2]))
    axdendro.set_xticks([])
    axdendro.set_yticks([])
    if condition == 0:
        pylab.title('Senescence', fontsize=20)
    elif condition == 1:
        pylab.title('Quiescence', fontsize=20)
    elif condition == 2:
        pylab.title('Proliferation', fontsize=20)
    else:
        raise ValueError
    
    #plot heatbar
    axheatbar = fig.add_axes([0.2,0.1,0.05,0.75])
    if condition < 3:
        label_color = Z['color_list']
        label_color.append(label_color[-1])
        label_color = np.array(label_color)
    axheatbar.set_xlim(0, 1)
    axheatbar.set_ylim(0, network_total)
    axheatbar.set_xticks([])
    axheatbar.set_yticks([])
    x = [0, 1]
    y = [[1, 1]] * network_total
    axheatbar.stackplot(x,y,colors=label_color)
    
    #plot distance matrix
    axmatrix = fig.add_axes([0.26,0.1,0.64,0.75])
    index = Z['leaves']
    attlistsorted = attlist[index,:]
    _, attlistlen = attlist.shape
    pchlder = int(attlistlen/3)
    if condition == 0:
        attlistsorted = attlistsorted[:,:pchlder]
    elif condition == 1:
        attlistsorted = attlistsorted[:,pchlder:2*pchlder]
    elif condition == 2:
        attlistsorted = attlistsorted[:,2*pchlder:]
    else:
        raise ValueError
    im = axmatrix.matshow(attlistsorted, aspect='auto', origin='lower', cmap=ListedColormap(['k', 'w']))
    axmatrix.set_yticks([])
    axmatrix.set_xticks(np.arange(len(node_list)))
    axmatrix.set_xticklabels(node_list, rotation=90)
    
    #display figure
    fig.show()
    
    if save == True:
        if condition == 0:
            shorthand = 'sene'
        elif condition == 1:
            shorthand = 'qui'
        elif condition == 2:
            shorthand = 'pro'
        else:
            raise ValueError
        if mark_pdk1:
            pylab.savefig('att_cluster_{}_{}.png'.format(shorthand, 'pdk1'))
        else:
            pylab.savefig('att_cluster_{}.png'.format(shorthand))
    
    if condition == 2:
        #mark PDK1 controlled models
        netparamslistsorted = np.array(netparamslist)[index]
        targetlistsorted = np.array(targetlist)[index]
        attractorlistsorted = np.array(attractorlist)[index]
        attractorlist_druggedsorted = np.array(attractorlist_drugged)[index]
        pdk1attractors_bf = []
        pdk1attractors_af = []
        pdk1_idxs = []#cluster position that needs to be yellow
        for idx, targets in enumerate(targetlistsorted):
            if ('PDK1',) in targets:
                pdk1_idxs.append(idx)
                pdk1attractors_bf.append(attractorlistsorted[idx])
                pdk1attractors_af.append(attractorlist_druggedsorted[idx][targets.index(('PDK1',))])
        pdk1_idxs = np.array(pdk1_idxs)
        with open('pdk1_attractors_bf.bin', 'wb') as f:
            pickle.dump(pdk1attractors_bf, f)
        with open('pdk1_attractors_af.bin', 'wb') as f:
            pickle.dump(pdk1attractors_af, f)
        
        #color_dict = {'b': '#0000FF', 'c': '#00EEEE', 'g': '#008000', 'k': , 'm': , 'r': '#FF0000', 'y': '#FFFF00'}
        label_color[pdk1_idxs] = 'y'
        link_colors = {}
        for i, leaves in enumerate(Y[:,:2].astype(int)):
            c1, c2 = (link_colors[leaf] if leaf > len(Y) else label_color[index.index(leaf)] for leaf in leaves)
            link_colors[i+1+len(Y)] = c1
    
    if condition >= 3:
        #plot only the pdk1 controllable models
        
        #plot dendrogram
        fig = pylab.figure()
        axdendro = fig.add_axes([0.09,0.1,0.1,0.75])
        if mark_pdk1:
            Z = sch.dendrogram(Y[pdk1_idxs,pdk1_idxs], orientation='left', color_threshold=None, link_color_func=lambda x: link_colors[x])
        else:
            Z = sch.dendrogram(Y[pdk1_idxs,pdk1_idxs], orientation='left', color_threshold=.5*max(Y[:,2]))
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        if condition == 0:
            pylab.title('Senescence', fontsize=20)
        elif condition == 1:
            pylab.title('Quiescence', fontsize=20)
        elif condition == 2:
            pylab.title('Proliferation', fontsize=20)
        else:
            raise ValueError
        
        #plot distance matrix
        axmatrix = fig.add_axes([0.2,0.1,0.7,0.75])
        index = Z['leaves']
        attlistsorted = attlist[index,:]
        attlistpdk1 = attlistsorted[pdk1_idxs,:]
        _, attlistlen = attlist.shape
        pchlder = int(attlistlen/3)
        if condition == 0:
            attlistpdk1 = attlistpdk1[:,:pchlder]
        elif condition == 1:
            attlistpdk1 = attlistpdk1[:,pchlder:2*pchlder]
        elif condition == 2:
            attlistpdk1 = attlistpdk1[:,2*pchlder:]
        else:
            raise ValueError
        im = axmatrix.matshow(attlistpdk1, aspect='auto', origin='lower', cmap=ListedColormap(['k', 'w']))
        axmatrix.set_yticks([])
        axmatrix.set_xticks(np.arange(len(node_list)))
        axmatrix.set_xticklabels(node_list, rotation=90)
        
        #display figure
        fig.show()
    
        if save == True:
            if condition == 0:
                shorthand = 'sene'
            elif condition == 1:
                shorthand = 'qui'
            elif condition == 2:
                shorthand = 'pro'
            else:
                raise ValueError
            if mark_pdk1:
                pylab.savefig('att_cluster_{}_{}.png'.format(shorthand, 'pdk1'))
            else:
                pylab.savefig('att_cluster_{}.png'.format(shorthand))
    
    if plot_tsne == True:
        data_ = attlistsorted
        tsne = TSNE(n_components=tsne_dim, random_state=0).fit_transform(data_)
        
        if tsne_dim == 2:
            plt.figure()
            plt.scatter(tsne[:, 0], tsne[:, 1], c=label_color, alpha=.2)
        elif tsne_dim == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(tsne[:, 0], tsne[:, 1], tsne[:, 2], c=label_color)
        else:
            raise ValueError
        
        if save == True:
            if mark_pdk1:
                pylab.savefig('att_cluster_tsne_{}_{}.png'.format(shorthand, 'pdk1'))
            else:
                pylab.savefig('att_cluster_tsne_{}.png'.format(shorthand))






















