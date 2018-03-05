import numpy as np
import pandas as pd
from PyGMO.problem import base

class OptWeightandBasal(base):
    """
    This is genetic algorithm for fitting weight and basal activity of a Boolean network.
    The algorithm takes network information in .sif and input and objective network conditions in a .csv format.
    The network information and objective condition is inherited as list and pd.DataFrame, respectively.
    The files with this code should provide an example for the format of the .csv file.
    The dimension of the parameter must be equal to the number of links(weights) + nodes(basal activities) - forced parameters.
    The parameters optimized by the genetic algorithm is a vector with link weights and node basal activities.
    """
	
    ###initial parameter setting
    def __init__(self):
        self._initialize()
        super(OptWeightandBasal, self).__init__(self.dim, self.i_dim)
        #super(OptWeightandBasal, self).__init__(self.dim, self.i_dim, self.num_conds) # for multivariate algorithms

        lb = np.zeros((self.dim,), dtype=np.double)
        ub = np.zeros((self.dim,), dtype=np.double)

        # Boundary of edges
        lb[:self.num_edges] = -self.wr #1 #-self.wr
        ub[:self.num_edges] = self.wr
        
        # Boundary of basal activity
        lb[self.num_edges:] = -self.wr
        ub[self.num_edges:] = self.wr

        self.set_bounds(lb, ub)

    def _initialize(self):

#        with open('netnodump.txt', 'r+') as netnodump:
#            netno = netnodump.read().strip()
#            netnodump.seek(0)
#            netnodump.write(str(int(netno) + 1))
#            netnodump.truncate()

        # netfilename = 'PKN_20_' + netno.strip() + '.sif'
        netfilename = "pkn_20_extra.sif"
        netdata = [] #list with network information
        netfile = open(netfilename, 'r')
        for link in netfile:
            netdata.append(tuple(link.strip().split('\t')))
        conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
        iptidx = []#indicies of the input conditions
        for idx, state in enumerate(list(conddata.iloc[0])):
            if state != 9:
                iptidx.append(idx)
        
        self.netdata = netdata
        self.conddata = conddata
        self.iptidx = iptidx
        self.netfilename = netfilename
        self.node_list = list(conddata.columns)
        self.num_conds = len(conddata.index) / 2 #number of objective conditions
        self.num_nodes = len(self.node_list)
        self.num_edges = len(self.netdata)
        
        #calculating maximum indegree to set the bounds below
        indeglist = {}
        for node in self.node_list:
            indeglist[node] = 0
        for link in self.netdata:
            indeglist[link[2]] += 1
        self.maxindeg = indeglist[max(indeglist, key = indeglist.get)]
        self.wr = 10 * (self.maxindeg + 2)
        
        #manually set parameters to force
        forcew = [['ULK1', 'ATG13', 2],
                  ['ULK1', 'ATG5', 2],
                  ['FOXO3', 'ATM/ATR', 1],
                  ['DNAdamage', 'ATM/ATR', 3],
                  ['DNAdamage_high', 'ATM/ATR_2', 2],
                  ['RB1', 'E2F1', 3],
                  ['ATM/ATR', 'E2F1', 1],
                  ['MAPK1', 'EIF4EBP1', 0],
                  ['MTOR', 'EIF4EBP1', 2],
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
                  ['IGF1R', 'KRAS', 1],#
                  ['ATP', 'ADP/ATP', 2],
                  ['lowNutrition', 'ADP/ATP', 1],
                  ['ATP', 'AMP/ATP', 2],
                  ['lowNutrition', 'AMP/ATP', 1],
                  ['TP53', 'AMPK', 1], 
                  ['ADP/ATP', 'AMPK', 1],
                  ['AMP/ATP', 'AMPK', 1],
                  ['lowNutrition', 'AMPK', 1],
                  ['Glycolysis', 'ATP', 1],
                  ['TCA cycle', 'ATP', 1],
                  ['BAX', 'CASP3', 1],
                  ['CASP9', 'CASP3', 2],
                  ['CDKN1A', 'CDK4', 2],
                  ['CDKN2A', 'CDK4', 2],
                  ['SIRT1', 'NFKB1', 1],
                  ['NFKBIE', 'NFKB1', 3],
                  ['NFKB1', 'NFKBIE', 1],
                  ['IKBKB', 'NFKBIE', 3],
                  ['MDM2', 'RB1', 1],
                  ['CDK4', 'RB1', 3],
                  ['ROS', 'DNAdamage', 0],
                  ['TP53', 'IGF1', 0],
                  ['INSR', 'IRS1', self.wr / 2],#basal activity set to 5to give half chance of influencing X; basal activity bound (1~9)
                  ['IGF1R', 'IRS1', 2],
                  ['S6K1', 'IRS1', 1],
                  ['KRAS', 'PIK3CA', 1],
                  ['IRS1', 'PIK3CA', 2]
                  ]#index and value of weight I want to manually adjust;['node1', 'node2', 1]
        forceb = [['ATG13', -1],
                  ['ATG5', -1],
                  ['ATM/ATR', -2],
                  ['ATM/ATR_2', -1],
                  ['E2F1', 1],
                  ['EIF4EBP1', 1],
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
                  ['KRAS', -1],#
                  ['ADP/ATP', 1],
                  ['AMP/ATP', 1],
                  ['AMPK', -2],
                  ['ATP', -1],
                  ['CASP3', -2],
                  ['CDK4', 1],
                  ['NFKB1', 2],
                  ['NFKBIE', 1],
                  ['RB1', 2],
                  ['IRS1', -1],
                  ['PIK3CA', -2]
                  ]#index and value of basal activity I want to manually adjust;['node1', 1]
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
        ##########initialize to no forced variables
        forcew = []
        forceb = []
        ##########temporary gesture
        dim = len(netdata) + len(conddata.columns) - len(forcew) - len(forceb)#genetic algorithm parameter dimension; number of edges(weights) + number of nodes(basal activities) - number of fixed parameters

        self.dim = dim
        self.i_dim = 0 #self.dim
        self.forcew = forcew
        self.forceb = forceb
        self.weight_matrix = np.zeros((self.num_nodes,
                                       self.num_nodes),
                                      dtype=np.int32)

        self.signs = np.ones((self.num_edges,), dtype=np.int32)

        self.ind_rows = np.zeros((self.num_edges,), dtype=np.int32)
        self.ind_cols = np.zeros((self.num_edges,), dtype=np.int32)
        
        for idx, edge in enumerate(self.netdata):
            if edge[1] == 'inhibit':
                self.signs[idx] = -1
            self.ind_rows[idx] = self.node_list.index(edge[0])
            self.ind_cols[idx] = self.node_list.index(edge[2])
        # end of for
        self.weight_matrix = np.zeros((self.num_nodes, self.num_nodes))

    def _index_hist(self, st, hist):
        for i, state in enumerate(hist):
            if (state == st).all():
                return i
        return -1
        
    ###the Boolean network model
    def _func(self, w, b, inistates):
        #self.weight_matrix[self.ind_rows, self.ind_cols] = self.signs*np.round(w)
        self.weight_matrix[self.ind_rows, self.ind_cols] = np.round(w)
        #self.weight_matrix = self.weight_matrix.astype(np.int64)
        #get the attractors
        state_history = [] # Record of state transition
        
        x_t1 = np.array(inistates)
        state_history.append(x_t1)
        while True:#state transition from t to t + 1
            
            x_t2 = np.dot(x_t1, self.weight_matrix) + b #transitioned state
            
            # Forcing the initial states
            for idx in self.iptidx:             
                if inistates[idx] == 0:
                    x_t2[idx] = 0
                elif inistates[idx] == 1:
                    x_t2[idx] = 1
            # end of for

            # Activation function
            x_t2[x_t2 > 0] = 1
            x_t2[x_t2 <= 0] = 0
            
            idx = self._index_hist(x_t2, state_history)
            if idx >= 0:
                # point attractor has list length of 1
                attractor = state_history[idx:]
                return attractor
            
            else: # move on
                state_history.append(x_t2.copy())
                x_t1 = x_t2
                

    ###reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):#params is edge weight and node basal level in that order
        w = np.array(params[:self.num_edges - len(self.forcew)])#weight parameters
        b = np.array(params[self.num_edges - len(self.forcew):])#basal activity parameter
		# inserting forced parameters
        for fw in self.forcew:#fw[0] is the index and fw[1] is the value to append
            w = np.insert(w, fw[0], fw[1])
        for fb in self.forceb:#fb[0] is the index and fb[1] is the value to append
            b = np.insert(b, fb[0], fb[1])
		
        fitness = 0.0
        for condno in range(self.num_conds):#for all of the conditions presented calculate the cost function and add all
            condno += 1
            inistates = self.conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
            #for idx in self.iptidx:#fixing input nodes using extreme basal activity settings
            #    if inistates[idx] == 0:
            #        b[idx] = -self.wr
            #    elif inistates[idx] == 1:
            #        b[idx] = self.wr #unnecessary redundancy
			
            inistates = np.array(inistates)
            inistates[inistates == 9] = 0
            predstateslist = self._func(w, b, inistates) #calculated states form the Boolean model
            objstates = np.array(self.conddata.loc['OUTPUT: ' + str(condno)].tolist())#objective states
            objboolmask = objstates != 9
            objphen = objstates[objboolmask]
            tmpfitlist = np.zeros((len(predstateslist),))
            for i, predstates in enumerate(predstateslist):
                predstates = np.array(predstates)
                predphen = predstates[objboolmask]
                tmpfit = np.sum((objphen - predphen)**2)
                tmpfitlist[i] = tmpfit
            fitness += np.mean(tmpfitlist)
        # end of for
        return (fitness,)

    #add an output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight and Basal activity Optimization for Boolean Network Model"
