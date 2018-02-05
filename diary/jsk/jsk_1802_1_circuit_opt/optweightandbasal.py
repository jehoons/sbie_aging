import numpy as np
import pandas as pd
from PyGMO.problem import base

#import network and objective network conditions with inputs
#import needs to be made here because of the problem with the pygmo package
#the number of iterations in the forbash.sh file needs to be consistent with the number of PKN_20 variants
with open('netnodump.txt', 'r+') as netnodump:
    netno = netnodump.read()
    netnodump.seek(0)
    netnodump.write(str(int(netno) + 1))
    netnodump.truncate
netfilename = 'PKN_20_' + netno + '.sif'
netdata = []#list with network information
netfile = open(netfilename, 'r')
for link in netfile:
    netdata.append(tuple(link.strip().split('\t')))
conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
iptidx = []#indicies of the input conditions
for idx, state in enumerate(list(conddata.iloc[0])):
    if state != 9:
        iptidx.append(idx)
forcew = [['ULK1', 'ATG13', 2],
          ['ULK1', 'ATG5', 2],
          ['FOXO3', 'ATM/ATR', 1],
          ['DNAdamage', 'ATM/ATR', 3],
          ['DNAdamage_high', 'ATM/ATR_2', 2],
          ['RB1', 'E2F1', 3],
          ['ATM/ATR', 'E2F1', 1],
          ['MAPK1', 'EIF4EBP1', 1],
          ['MTOR', 'EIF4EBP1', 3],
          ['FOXO3', 'G6PC', 2],
          ['SLC2A4', 'Glucose', 2],
          ['ATM/ATR_2', 'HIPK2', 2],
          ['IGF1', 'IGF1R', 3],
          ['TP53', 'IGF1R', 1],
          ['TP53_s46', 'IGF1R', 0],
          ['NFKB1', 'IL1B', 2],
          ['NFKB1', 'IL6', 2],
          ['lowNutrition', 'INSR', 2],
          ['MAP2K3/MAP2K6', 'MAPK14', 2],
          ['AMPK', 'NAD+', 2],
          ['TSC2', 'RHEB', 2],
          ['MTOR', 'S6K1', 2],
          ['FOXO3', 'SOD2', 2],
          ['Glycolysis', 'TCA cycle', 2],
          ['NFKB1', 'TNF', 2],
          ['AMPK', 'ULK1', 1],
          ['MTOR', 'ULK1', 3],
          ['INSR', 'KRAS', 1],
          ['IGF1R', 'KRAS', 1]]#index and value of weight I want to manually adjust;['node1', 'node2', 1]
forceb = [['ATG13', -1],
          ['ATG5', -1],
          ['ATM/ATR', -2],
          ['ATM/ATR_2', -1],
          ['E2F1', 1],
          ['EIF4EBP1', 2],
          ['G6PC', -1],
          ['Glucose', -1],
          ['HIPK2', -1],
          ['IGF1R', -1],
          ['IL1B', -1],
          ['IL6', -1],
          ['INSR', 1],
          ['MAPK14', -1],
          ['NAD+', -1],
          ['RHEB', 1],
          ['S6K1', -1],
          ['SOD2', -1],
          ['TCA cycle', -1],
          ['TNF', -1],
          ['ULK1', 2],
          ['KRAS', -1]]#index and value of basal activity I want to manually adjust;['node1', 1]
tmpforcew = []#set appropriate format for forcew
for fw in forcew:
    if (fw[0], 'activate', fw[1]) in netdata:
        tmpforcew.append((netdata.index((fw[0], 'activate', fw[1])), fw[2]))
    elif ((fw[0], 'inhibit', fw[1]) in netdata):
        tmpforcew.append((netdata.index((fw[0], 'inhibit', fw[1])), fw[2]))
forcew = sorted(tmpforcew)
forceb = [(list(conddata.columns).index(fb[0]), fb[1]) for fb in forceb]#set appropriate format for forceb
for idx in iptidx:#integrating input and force basal activities together for technical reasons; both are essentially basal value fixation
    forceb.append((idx, 9))
forceb = sorted(forceb)
dim = len(netdata) + len(conddata.columns) - len(forcew) - len(forceb)#genetic algorithm parameter dimension; number of edges(weights) + number of nodes(basal activities) - number of fixed parameters

class OptWeightandBasal(base):
    """
    This is genetic algorithm for fitting weight and basal activity of a Boolean network.
    The algorithm takes network information in .sif and input and objective network conditions in a .csv format.
    The network information and objective condition is inherited as list and pd.DataFrame, respectively.
    The files with this code should provide an example for the format of the .csv file.
    The dimension of the parameter must be equal to the number of links(weights) + nodes(basal activities).
    The parameters optimized by the genetic algorithm is a vector with link weights and node basal activities in that order
    Because the weights and basal activities are optimized on a single vector,
    manual adjustments have been made for basal activity boundaries; they have been repositioned by maxindegree
    """
    
    ###initial parameter setting
    def __init__(self, dim = dim, i_dim = dim, netdata = netdata, conddata = conddata, iptidx = iptidx, forcew = forcew, forceb = forceb, netfilename = netfilename):
        super(OptWeightandBasal, self).__init__(dim, i_dim)
        self.netdata = netdata
        self.conddata = conddata
        self.iptidx = iptidx
        self.netfilename = netfilename
        self.node_list = list(conddata.columns)
        self.num_conds = len(conddata.index) / 2#number of objective conditions
        self.num_nodes = len(self.node_list)
        self.num_edges = len(self.netdata)
        self.forcew = forcew
        self.forceb = forceb
        
        #calculating maximum indegree to set the bounds below
        indeglist = {}
        for node in self.node_list:
            indeglist[node] = 0
        for link in self.netdata:
            indeglist[link[2]] += 1
        self.maxindeg = indeglist[max(indeglist, key = indeglist.get)]

        #set the bounds
        #the boundaries must encompass the maximum indegree and minimum basal level
        #as the basal level in included with the edge weight in the parameter,
        #the basal level's bounds has to be adjusted within the _func function by maxindegree
        lb = np.zeros((dim,), dtype=np.double)
        ub = np.zeros((dim,), dtype=np.double)
        self.wr = 2 * (self.maxindeg + 1 + 1)#weight resolution; maxindegree + basal activity + difference
        lb[:] = 1
        ub[:] = self.wr
        self.set_bounds(lb, ub)#the double is for manual tweak to cover the variation for the basal acitvities

    ###the Boolean network model
    def _func(self, w, b, inistates):
        #generate the weight matrix
        weightmat = np.zeros((self.num_nodes, self.num_nodes), dtype = np.int32)
        for idx, edge in enumerate(self.netdata):
            if edge[1] == 'inhibit':
                w[idx] = -w[idx]
            row = self.node_list.index(edge[0])
            col = self.node_list.index(edge[2])
            weightmat[row][col] = w[idx]
        weightmat = np.vstack((weightmat, b))#stack basal activity levels
        self.weight_matrix = weightmat
        #get the attractors
        state_history = []#record of state transition
        while 1:#state transition from t to t + 1
            initial_state = list(inistates)#save initial states of nodes at t
            inistates.append(1)#attaching basal activity levels at the end; complementary with stacked basal activity on the weight matrix
            inistates = np.array(inistates, dtype = np.int64)
            self.weight_matrix.astype(np.int64)
            transtates = np.dot(inistates, self.weight_matrix)#transitioned state
            transtates[transtates > 0] = 1
            transtates[transtates <= 0] = 0
            transed_state = list(transtates)#temporarily save states of nodes at t + 1
            state_history.append(initial_state)#save the transient state

            if transed_state in state_history:#get the point attractor or the cyclic attractor; list of the state trajectory since the redundant state
                self.state_history = state_history
                attractor = state_history[state_history.index(transed_state):]#point attractor has list length of 1
                return attractor
            
            else:#move on
                inistates = list(transed_state)

    ###reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):#params is edge weight and node basal level in that order
        w = np.array(params[:self.num_edges - len(self.forcew)])#weight parameters
        w = ((w + 1) / 2).astype(np.int32)#adjusting the boundaries to suit the weight parameters
        b = np.array(params[self.num_edges - len(self.forcew):])#basal activity parameter
        b = np.array([a - self.wr / 2 if a > self.wr / 2 else a - self.wr / 2 - 1 for a in b]).astype(np.int32)#adjusting the boundaries to cover negative basal activity levels
        for fw in self.forcew:#fw[0] is the index and fw[1] is the value to append
            w = np.insert(w, fw[0], fw[1])
        for fb in self.forceb:#fb[0] is the index and fb[1] is the value to append
            b = np.insert(b, fb[0], fb[1])
        fitness = 0.0
        for condno in range(self.num_conds):#for all of the conditions presented calculate the cost function and add all
            condno += 1
            inistates = self.conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
            for idx in self.iptidx:#fixing input nodes using extreme basal activity settings
                if inistates[idx] == 0:
                    b[idx] = -self.maxindeg
                elif inistates[idx] == 1:
                    b[idx] = self.maxindeg
            inistates = list(map(lambda x: int(str(x).replace('9', '0')), inistates))
            predstateslist = self._func(w, b, inistates)#calculated states form the Boolean model
            objstates = np.array(self.conddata.loc['OUTPUT: ' + str(condno)].tolist())#objective states
            objboolmask = np.array([state != 9 for state in objstates])#Boolean mask for selecting only the states required for cost function calculation
            objphen = list(objstates[objboolmask])
            tmpfitlist = []#list of fitness value for all of the states within a cyclic attractor; only one of if a point attractor
            for predstates in predstateslist:
                predstates = np.array(predstates)
                predphen = list(predstates[objboolmask])
                tmpfit = sum((objphen[idx] - predphen[idx])**2 for idx in range(len(objphen)))#fitness function calulated for this condition and added
                tmpfitlist.append(tmpfit)
            fitness += np.mean(np.array(tmpfitlist))#fitness is the mean of fitness of all the states in a cyclic attractor
        return (fitness,)

    #add an output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight and Basal activity Optimization for Boolean Network Model"
