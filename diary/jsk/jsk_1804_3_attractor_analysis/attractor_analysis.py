import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from matplotlib.colors import ListedColormap

#import
netfilename = "PKN_22.sif"
netdata = [] #list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

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

paramsdumpname = "netparamsdump.txt"
paramsdata = []
paramsfile = open(paramsdumpname, 'r')
for params in paramsfile:
    paramsdata.append(params.strip().split('\t'))

totdatano = int((len(paramsdata) - 4) / 3)
ws = np.zeros((totdatano, len(netdata)), dtype=np.float64)
bs = np.zeros((totdatano, len(node_list) - len(iptidx)), dtype=np.float64)
params = []

for datano in range(totdatano):
    ws[datano] = np.array(paramsdata[datano * 3 + 4], dtype=np.float64)
    bs[datano] = np.array(paramsdata[datano * 3 + 5], dtype=np.float64)
    params.append(np.concatenate((ws[datano], bs[datano])))
params = np.array(params)

tempbs = np.transpose(bs).copy()

bs = np.zeros((len(node_list), totdatano), dtype=np.float64)
i = 0
for idx in range(len(node_list)):
    if idx in iptidx:
        bs[idx] = np.array([99] * totdatano)
    else:
        bs[idx] = tempbs[i]
        i += 1
bs = np.transpose(bs)

netparamslist = []
for w, b in zip(ws, bs):
    netparamslist.append((w, b))

bsdf = pd.DataFrame(bs)
wsdf = pd.DataFrame(ws)


#simulator
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


def _index_hist(st, hist):
    for i, state in enumerate(hist):
        if (state == st).all():
            return i
    return -1
    
###the Boolean network model
def _func(w, b, inistates):
    weight_matrix[ind_rows, ind_cols] = signs*np.round(w)
    #self.weight_matrix = self.weight_matrix.astype(np.int64)
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
attlist = np.array(attlist)


###############################################################################
###clustering analysis
import pylab
import scipy.cluster.hierarchy as sch

data = pd.DataFrame(attlist)
data = data.T

#generate correlation matrix
corrmat = pd.DataFrame(0, index = data.columns, columns = data.columns).astype('float64')
for gene1 in data.columns:
    for gene2 in data.columns:
        #corrmat[gene1][gene2] = data[gene1].corr(data[gene2])
        corrmat[gene1][gene2] = sum(data[gene1] - data[gene2])

# Compute and plot dendrogram
D=corrmat.values
D=np.nan_to_num(D)
fig = pylab.figure()
axdendro = fig.add_axes([0.09,0.1,0.1,0.75])
Y = sch.linkage(D, method='centroid')
Z = sch.dendrogram(Y, orientation='left')
axdendro.set_xticks([])
axdendro.set_yticks([])
#pylab.title('Senescence', fontsize=20)
#pylab.title('Quiescence', fontsize=20)
pylab.title('Proliferation', fontsize=20)

# Plot distance matrix.
axmatrix = fig.add_axes([0.2,0.1,0.7,0.75])
index = Z['leaves']
#D = D[index,:]
#D = D[:,index]
#im = axmatrix.matshow(D, aspect='auto', origin='lower')
attlistsorted = attlist[index,:]
#attlistsorted = attlistsorted[:,:42]
#attlistsorted = attlistsorted[:,42:84]
attlistsorted = attlistsorted[:,84:]
im = axmatrix.matshow(attlistsorted, aspect='auto', origin='lower', cmap=ListedColormap(['k', 'w']))
axmatrix.set_yticks([])
axmatrix.set_xticks(np.arange(len(node_list)))
axmatrix.set_xticklabels(node_list, rotation=90)


# Display figure.
fig.show()
#pylab.savefig('corr.png')

#pd.DataFrame(attlistsorted[:,:42], columns=node_list).to_csv('attlist_sene.csv')



'''
clusterfilename = "tcluster1.txt"
clustlabel1= []
clusterfile = open(clusterfilename, 'r')
for shit in clusterfile:
    clustlabel1.append(shit.strip().split())
clustlabel1 = clustlabel1[1:]
clusterfilename = "tcluster2.txt"
clustlabel2= []
clusterfile = open(clusterfilename, 'r')
for shit in clusterfile:
    clustlabel2.append(shit.strip().split())
clustlabel2 = clustlabel2[1:]
'''

'''
#actual clustering
num = 2
data = attlist
kmeans = KMeans(n_clusters=num, random_state=0).fit_predict(data)
tsne = TSNE(n_components=num, random_state=0).fit_transform(data)


#LABEL_COLOR_MAP = {'-1' : '0', '1' : 'tab:blue', '2' : 'tab:orange', '3' : 'tab:green',
#                   '4' : 'tab:red', '5' : 'tab:purple', '6' : 'tab:brown',
#                   '7' : 'tab:pink', '8' : 'tab:gray', '9' : 'tab:olive',
#                   '10' : 'tab:cyan'}
#label_color = [LABEL_COLOR_MAP[l[1]] for l in clustlabel2]
LABEL_COLOR_MAP = {0 : 'r',1 : 'b',2 : 'y'}
label_color = [LABEL_COLOR_MAP[l] for l in kmeans]
if num==2:
    plt.figure()
    plt.scatter(tsne[:, 0], tsne[:, 1], c=label_color, alpha=.5)
elif num==3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(tsne[:, 0], tsne[:, 1], tsne[:, 2], c=label_color)


###correlation clustering analysis
import pylab
import scipy.cluster.hierarchy as sch

def check(cno, tmpcluster):
    if cno < geneno:
        tmpcluster.append(int(cno))
    else:
        Yidx = int(cno - geneno)
        tmpcluster = check(Y[Yidx][0], tmpcluster)
        tmpcluster = check(Y[Yidx][1], tmpcluster)
    return tmpcluster

params = pd.DataFrame(params)

#generate correlation matrix
corrmat = pd.DataFrame(0, index = params.columns, columns = params.columns).astype('float64')
for gene1 in params.columns:
    for gene2 in params.columns:
        corrmat[gene1][gene2] = params[gene1].corr(params[gene2])

# Compute and plot dendrogram
D=corrmat.values
fig = pylab.figure()
axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
Y = sch.linkage(D, method='centroid')
Z = sch.dendrogram(Y, orientation='left')
axdendro.set_xticks([])
axdendro.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
index = Z['leaves']
D = D[index,:]
D = D[:,index]
im = axmatrix.matshow(D, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
pylab.colorbar(im, cax=axcolor)

# Display figure.
fig.show()
#pylab.savefig('corr.png')

#despense clusters
geneno = len(params.columns)
clusterno = 2#number of clusters
clusterYno = []
for idx in range(1, clusterno):
    clusterYno.append(Y[-idx][0])
    clusterYno.append(Y[-idx][1])
for i in range(clusterno - 2):
    clusterYno.remove(max(clusterYno))
clusters = []
for cno in clusterYno:
    tmpcluster = []
    tmpcluster = check(cno, tmpcluster)
    clusters.append(list(np.array(params.columns)[tmpcluster]))


stuff=['_s', '_q', '_p']
colname = []
for st in stuff:
    for node in node_list:
        colname.append(node + st)
#pd.DataFrame(np.unique(attlist, axis=0), columns=colname).to_csv('attlist.csv')
'''









