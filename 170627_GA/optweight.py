import numpy as np
from PyGMO.problem import base
import numpy as np

class OptWeight(base):
    """
    This GA is for weight fitting for boolean network.
    You have to have node, edge, basal level and attractor information.
    See the example files ('node_info.txt','edge_info.txt','obj_attractor.txt', ??basal)
    and prepare your own files in same format.
    Dimesion must be equal to the number of links.
    """
    ### initial parameter setting
    def __init__(self, dim = 24, i_dim=24):  # num of link 13(weight) + 11(basal)
        super(OptWeight, self).__init__(dim, i_dim)

        # attractor Data
        self.data = []
        f = open('obj_attractor.txt', 'r')
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
        self.basal = [0] * self.n_node
        for line in file_edge:
            self.edge_list.append(line.strip().split(' '))
        self.n_weight = len(self.edge_list)
        try:
            len(self.edge_list)+len(self.basal) == dim
        except ValueError:
            print("Oops! Dimension was wrong. Dimension must be equal to the num of link")
        else:
            pass

        # set the bounds
        lb = np.zeros((dim,), dtype=np.double)
        ub = np.zeros((dim,), dtype=np.double)
        lb[:] = -10
        for i,j in enumerate(self.edge_list):
            if j[1] == '((activate))':
                lb[i] = 0
        ub[:] = 10  # ub= 1 + max(degree) & larger than basal
        for i,j in enumerate(self.edge_list):
            if j[1] == '((inhibit))':
                ub[i] = 0
        self.set_bounds(lb, ub)


    ### Model
    def _func(self, w, k):
        # make weight matrix
        a = np.zeros((self.n_node, self.n_node))
        y = 0
        for (enum, e) in enumerate(self.edge_list):
            row = self.node_list.index(e[0])
            col = self.node_list.index(e[2])
            a[row][col] = w[y]
            y += 1
        a = np.vstack((a, w[-11:]))
        self.weight_matrix = a

        # get attractor
        dic_state = {}
        while 1:
            # transition from t to t+1
            x_sum = k
            initial_state = tuple(k)  # state of 3 nodes at initial t
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
    def _objfun_impl(self, params):
        w = np.round(params)
        s = 0.0
        for state in self.data:
            initial = list(state[0])
            k = [int(i) for i in initial]
            # ex) k= [1 1 0 1 1]
		
            obj = list(state[1])
            obj_state = [int(j) for j in obj]
            # ex) obj_state = [1, 1, 1, 1, 1]

            y = self._func(w, k)
            # ex) y = [1, 0, 1, 1, 1]
            advantage = 0
            for i in w:
                if int(i) == 0:
                    advantage+=1
            constant = 5
            s += sum((obj_state[i] - y[i])**2 for i in range(len(y)))+constant*advantage

        return (s,)

    # Add some output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight optimization for boolean toy model"