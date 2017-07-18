import numpy as np
from PyGMO.problem import base
import numpy as np
import json

class OptWeight(base):
    """
    This GA is for weight fitting for boolean network.
    You have to have node, edge, basal level and attractor information.
    See the example files ('node_info.txt','edge_info.txt','obj_attractor.txt', ??basal)
    and prepare your own files in same format.
    Dimesion must be equal to the number of links.
    """
    ### initial parameter setting
    def __init__(self, dim = 10, i_dim=10):  # num of link 13(weight) + 11(basal)
        super(OptWeight, self).__init__(dim,i_dim)
        with open('data.txt', 'r') as outfile:
            J = json.load(outfile)
        # simulation condition
        # self.sample_size = J['simulation_condition']['sample_size']
        self.repeat = J['simulation_condition']['repeat']
        self.lam = J['simulation_condition']['lam']
        self.step_limit = J['simulation_condition']['step_limit']
        # input condition
        self.input_dic = {}
        for input_condition in J['input_condition']:
            if input_condition:
                break
            self.input_dic[input_condition['node']]= input_condition['state']
        # node
        self.node_list = []
        for node in J['Node']:
            self.node_list.append(node['name'])
        self.n_node = len(self.node_list)
        basal = [0] * self.n_node
        # edge
        self.edge_list = []
        for pkn in J['PKN']:
            self.edge_list.append([pkn['source'],pkn['interaction'],pkn['target']])
        try:
            len(self.edge_list)+len(basal) == dim
        except ValueError:
            print("Oops! Dimension was wrong. Dimension must be equal to the num of parameters")
        else:
            pass
        # train graph Data
        self.train_data = []
        for train_graph in J['train_graph']:
            self.train_data.append([train_graph['before'], train_graph['pert'], train_graph['state'], train_graph['after']])
        # set the bounds
        lb = np.zeros((dim,), dtype=np.double)
        ub = np.zeros((dim,), dtype=np.double)
        lb[:] = -2
        for i,j in enumerate(self.edge_list):
            if j[1] == '((activate))':
                lb[i] = 0
        ub[:] = 2
        for i,j in enumerate(self.edge_list):
            if j[1] == '((inhibit))':
                ub[i] = 0
        self.set_bounds(lb, ub)

    def get_attractor(self, Weight, initial_state, state_traj, pert_node=-1, state = 0):
        # transition from t to t+1
        x = initial_state[:]
        if pert_node != -1:
            x[pert_node] = int(state)
        x.append(1)  # add basal level position
        x_array = np.array(x)
        traj_array = np.dot(Weight, x_array.transpose())

        traj_array[traj_array > 0] = int(1)
        traj_array[traj_array <= 0] = int(0)
        traj_list = [int(x) for x in traj_array]
        if pert_node != -1:
            traj_list[pert_node] = int(state)
        state_traj.append(initial_state)
        if traj_list == initial_state:
            return traj_list
        elif traj_list in state_traj:
            return state_traj[state_traj.index(traj_list):]
            # return traj_list
        else:
            return self.get_attractor(Weight=Weight, initial_state=traj_list, state_traj=state_traj, pert_node = pert_node, state=state)

    ### Model
    def _func(self, w, initial_attractor, pert=-1,state=0):
        # make weight matrix
        W_matrix = np.zeros((self.n_node, self.n_node))
        for (enum, e) in enumerate(self.edge_list):
            source = self.node_list.index(e[0])
            target = self.node_list.index(e[2])
            W_matrix[target][source] = w[enum]
        basal_array = np.array(w[-self.n_node:])
        basal_array.shape = (self.n_node, 1)
        W_matrix = np.hstack((W_matrix, basal_array))
        state_traj = []
        result = self.get_attractor(Weight=W_matrix, initial_state=list(initial_attractor), state_traj=state_traj, pert_node=pert, state=state)
        return result

    # Reimplement the virtual method that defines the objective function
    def _objfun_impl(self, params):
        w = list(np.round(params))
        total_s = 0.0
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
            result_att = self._func(w, initial)
            if type(result_att[0]) is int:
                attractor_set.append(result_att)
            elif type(result_att[0]) is list:
                for n, att in enumerate(result_att):
                    attractor_set.append(att)
        #make reachability graph
        reachability = []
        for initial_att in attractor_set:
            for i in range(self.n_node):
                final_att = self._func(w, initial_att, pert=int(i), state = 0)
                if type(final_att[0]) is int:
                    reachability.append([initial_att,i,0,final_att])
                elif type(final_att[0]) is list:
                    for att in final_att:
                        reachability.append([initial_att,i,0, att])
                final_att = self._func(w, initial_att, pert=int(i), state=1)
                if type(final_att[0]) is int:
                    reachability.append([initial_att,i,1,final_att])
                elif type(final_att[0]) is list:
                    for att in final_att:
                        reachability.append([initial_att,i,1, att])
        for state in self.train_data:
            obj = [[int(x) for x in list(state[0])],[int(x) for x in list(state[3])]]
            s = []
            for sim in reachability:
                if int(sim[1]) == int(state[1])&int(sim[2]) == int(state[2]):
                    s1 = [(int(i)-int(j))**2 for i, j in zip(obj[0],sim[0])]
                    s2 = [(int(i)-int(j))**2 for i, j in zip(obj[1],sim[3])]
                    s.append(sum(s1)+sum(s2))
            total_s += int(min(s))
        adv = 0
        for x in w:
            if int(x) == 0:
                adv+=1
        total_s = total_s - self.lam*adv
        return (total_s,)

    # Add some output to __repr__
    def human_readable_extra(self):
        return "\n\tWeight optimization for boolean toy model"