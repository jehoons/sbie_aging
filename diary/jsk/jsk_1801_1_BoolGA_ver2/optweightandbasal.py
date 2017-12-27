import numpy as np
import pandas as pd
from PyGMO.problem import base

#import network and objective network conditions with inputs as default
conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
netdata = []#list with network information
netfile = open('PKN_19.sif', 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))

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
    def __init__(self, dim = 1, i_dim = 1, netdata = netdata, conddata = conddata):
        super(OptWeightandBasal, self).__init__(dim, i_dim)
        self.netdata = netdata
        self.conddata = conddata
        self.node_list = list(conddata.columns)
        self.num_conds = len(conddata.index) / 2#number of objective conditions
        self.num_nodes = len(self.node_list)
        self.num_edges = len(self.netdata)
        
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
        lb[:] = 1
        ub[:] = 2 * self.maxindeg#weight resolution; set as maxindegree(single link should be possible to overthrow all the other)
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

            if transed_state in state_history:#get the point attractor and cyclic attractor; list of state trajectory since the redundant state
                self.state_history = state_history
                attractor = state_history[state_history.index(transed_state):]#point attractor has list length of 1
                return attractor
            
            else:#move on
                inistates = list(transed_state)

    ###reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):#params is edge weight and node basal level in that order
        w = np.array(params[:self.num_edges])#weight parameters
        w = ((w + 1) / 2).astype(np.int32)#adjusting the boundaries to suit the weight parameters
        b = np.array(params[self.num_edges:])#basal activity parameter
        bmask = [1 if s % 2 == 1 else -1 for s in b]
        b = ((b + 1) / 2).astype(np.int32)
        b = np.multiply(b, bmask)#adjusting the boundaries to cover negative basal activity levels
        fitness = 0.0
        
        print('new condition')
        
        for condno in range(self.num_conds):#for all of the conditions presented calculate the cost function and add all
            condno += 1
            inistates = self.conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
            predstateslist = self._func(w, b, inistates)#calculated states form the Boolean model
            
            print('new attractor')
            print(predstateslist)
            
            objstates = np.array(self.conddata.loc['OUTPUT: ' + str(condno)].tolist())#objective states
            objboolmask = np.array([state != 9 for state in objstates])#Boolean mask for selecting only the states required for cost function calculation
            objphen = list(objstates[objboolmask])
            tmpfitlist = []#list of fitness value for all of the states within a cyclic attractor; only one of if a point attractor
            for predstates in predstateslist:
                predstates = np.array(predstates)
                predphen = list(predstates[objboolmask])
                tmpfit = sum((objphen[idx] - predphen[idx])**2 for idx in range(len(objphen)))#fitness function calulated for this condition and added
                tmpfitlist.append(tmpfit)
            
            print('new phenotype')
            print(tmpfitlist)
            
            fitness += np.mean(np.array(tmpfitlist))#fitness is the mean of fitness of all the states in a cyclic attractor
            
            print('new fitness')
            print(fitness)
            
        return (fitness,)

    #add an output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight and Basal activity Optimization for Boolean Network Model"
