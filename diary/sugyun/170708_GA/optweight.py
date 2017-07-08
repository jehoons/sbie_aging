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
    def __init__(self, dim = 9, i_dim=9):  # num of link 13(weight) + 11(basal)
        super(OptWeight, self).__init__(dim, i_dim)

        # attractor Data
        self.data = []
        f = open('train_graph.txt', 'r')
        for n, line in enumerate(f):
            if n == 0:
                continue
            (ini, pert, obj) = line.strip().split('\t')
            self.data.append((ini, pert, obj))

        # edge & node & basal level Data
        file_edge = open('edge_info.txt', 'r')
        file_node = open('node_info.txt', 'r')
        self.node_list = []
        self.edge_list = []
        for line in file_node:
            self.node_list.append(line.strip())
        self.n_node = len(self.node_list)
        basal = [0] * self.n_node
        for line in file_edge:
            self.edge_list.append(line.strip().split(' '))
        try:
            len(self.edge_list)+len(basal) == dim
        except ValueError:
            print("Oops! Dimension was wrong. Dimension must be equal to the num of parameters")
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

    def get_attractor(self, Weight, initial_input, state_traj, pert_node=-1):
        # transition from t to t+1
        x = initial_input[:]
        if pert_node != -1:
            x[pert_node] = 0
        # x[1] = 0
        # x[2] = 0
        # x[3] = 0
        x.append(1)  # add basal level position
        x_array = np.array(x)
        traj_array = np.dot(Weight, x_array.transpose())

        traj_array[traj_array > 0] = int(1)
        traj_array[traj_array <= 0] = int(0)
        traj_list = [int(x) for x in traj_array]
        if pert_node != -1:
            traj_list[pert_node] = 0
        # traj_list[1] = 0
        # traj_list[2] = 0
        # traj_list[3] = 0
        state_traj.append(initial_input)
        if traj_list == initial_input:
            return traj_list
        elif traj_list in state_traj:
            return state_traj[state_traj.index(traj_list):]
            # return traj_list
        else:
            return self.get_attractor(Weight=Weight, initial_input=traj_list, state_traj=state_traj, pert_node = pert_node)

    ### Model
    def _func(self, w, k, pert=-1):
        # make weight matrix
        W = np.zeros((self.n_node, self.n_node))
        for (enum, e) in enumerate(self.edge_list):
            source = self.node_list.index(e[0])
            target = self.node_list.index(e[2])
            W[target][source] = w[enum]
        basal_array = np.array(w[-self.n_node:])
        basal_array.shape = (self.n_node, 1)
        W = np.hstack((W, basal_array))
        state_traj = []
        result = self.get_attractor(Weight=W, initial_input=list(k), state_traj=state_traj, pert_node=pert)
        return result

    # Reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):
        w = list(np.round(params))
        s = 0.0
        #all initial pts
        initial_num = 2 ** self.n_node
        initial_bin = []
        for i in range(initial_num):
            initial_bin.append(bin(i)[2:])
        initial_state = []
        for i in initial_bin:
            default = '0' * (self.n_node - len(i)) + i
            initial_state.append([int(s) for s in default])
        #get unique attractor
        attractor_set = []
        for initial in initial_state:
            y = self._func(w, initial)
            if type(y[0]) is int:
                attractor_set.append(y)
            elif type(y[0]) is list:
                for n, att in enumerate(y):
                    attractor_set.append(att)
        #make reachability graph
        reachability = []
        for k in attractor_set:
            for i in range(self.n_node):
                y = self._func(w, k, pert=int(i))
            # ex) y = [1, 0, 1, 1, 1]
                if type(y[0]) is int:
                    reachability.append([k,i,y])
                elif type(y[0]) is list:
                    for att in y:
                        reachability.append([k, i, att])
        total_s = 0
        for state in self.data:
            obj = [list(state[0]),list(state[2])]
            s = []
            for sim in reachability:
                if int(sim[1]) == int(state[1]):
                    s1 = [(int(i)-int(j))**2 for i, j in zip(obj[0],sim[0])]
                    s2 = [(int(i)-int(j))**2 for i, j in zip(obj[1],sim[2])]
                    s.append(sum(s1)+sum(s2))
            total_s += int(min(s))
        return (total_s,)

    # Add some output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight optimization for boolean toy model"