import numpy as np
from PyGMO.problem import base

class OptWeight(base):
    """
    This GA is for weight fitting for boolean network.
    You have to have node, edge, basal level and attractor information.
    See the example files ('node_info.txt','edge_info.txt','obj_attractor.txt', ??basal)
    and prepare your own files in same format.
    Dimesion must be equal to the number of links.
    """
    ### initial parameter setting
    def __init__(self, dim = 128, i_dim=128):#edgeno + nodeno
        super(OptWeight, self).__init__(dim, i_dim)

        # attractor Data
        self.data = []
        f = open('obj_attractor_sample.txt', 'r')###need to change this later
        for line in f:
            (ini, obj) = line.strip().split('\t')
            self.data.append((ini, obj))

        # edge & node & basal level Data
        file_edge = open('edge_info.txt', 'r')
        file_node = open('node_info.txt', 'r')
        self.node_list = []
        self.edge_list = []
        for line in file_node:
            self.node_list.append(line.strip())
        self.n_node = len(self.node_list)
        for line in file_edge:
            self.edge_list.append(line.strip().split(' '))
        try:
            len(self.edge_list) == dim
        except ValueError:
            print("Oops! Dimension was wrong. Dimension must be equal to the num of link")
        else:
            pass

        #set the bounds
        #the boundaries must encompass the maximum indegree and minimum basal level
        #as the basal level in included with the edge weight in the parameter,
        #the basal level's bounds has to be set manually within the _func function
        lb = np.zeros((dim,), dtype=np.double)
        ub = np.zeros((dim,), dtype=np.double)
        self.w_res = 5#max indegree
        lb[:] = 1
        ub[:] = self.w_res
        self.set_bounds(lb, ub)


    ### Model
    def _func(self, w, k, b):
        w = list(w)
        # make weight matrix
        a = np.zeros((self.n_node, self.n_node))
        y = 0
        for (enum, e) in enumerate(self.edge_list):
            if e[1] == '((inhibit))':
                if w[y] > 0:
                    w[y] = -w[y]
            row = self.node_list.index(e[0])
            col = self.node_list.index(e[2])
            a[row][col] = w[y]
            y += 1
        a = np.vstack((a, np.array(b)))#stack basal level
        self.weight_matrix = a

        # get attractor
        dic_state = {}
        while 1:
            # transition from t to t+1
            x_sum = k
            initial_state = tuple(k)  # state of nodes at initial t
            x_sum.append(1)  # add basal level
            x_sum = np.array(x_sum)
            traj_array = np.dot(x_sum, self.weight_matrix)
            traj_array[traj_array > 0] = 1
            traj_array[traj_array <= 0] = 0
            traj_array_tuple = tuple(int(i) for i in traj_array)

            # get fixed attractor
            if initial_state == traj_array_tuple:
                return list(traj_array_tuple)

            # get cyclic attractor
            elif traj_array_tuple in dic_state:
                return list(traj_array_tuple)
            dic_state[initial_state] = 0
            k = list(traj_array_tuple)

    ### Reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):#params is edge weight and node basal level in that order
        w = params[0:78]
        b = list(np.array(params[78:128]) - 2)#boundary repositioning for basal activites
        s = 0.0
        for state in self.data:
            initial = list(state[0])
            k = [int(i) for i in initial]
            y = self._func(w, k, b)
            y = np.array(y)
            
            obj = np.array(list(state[1]))
            obj_bool = np.array([i!='-' for i in obj])
            pred_state = list(y[obj_bool])
            obj_state = list(obj[obj_bool])
            obj_state = [int(i) for i in obj_state]
            
            s += sum((obj_state[i] - pred_state[i])**2 for i in range(len(obj_state)))
            
        return (s,)

    # Add some output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight optimization for boolean toy model"