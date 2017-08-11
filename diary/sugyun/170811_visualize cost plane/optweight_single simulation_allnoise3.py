import numpy as np
import json
with open('data_noise.txt', 'r') as outfile:
    J = json.load(outfile)
# simulation condition
sample_size = J['simulation_condition']['sample_size']
repeat = J['simulation_condition']['repeat']
lam = J['simulation_condition']['lam']
step_limit = J['simulation_condition']['step_limit']
# input condition
input_dic = {}
for input_condition in J['input_condition']:
    if input_condition:
        break
    input_dic[input_condition['node']] = input_condition['state']
# node
node_list = []
for node in J['Node']:
    node_list.append(node['name'])
n_node = len(node_list)
basal = [0] * n_node
# edge
edge_list = []
for pkn in J['PKN']:
    edge_list.append([pkn['source'], pkn['interaction'], pkn['target']])

# train graph Data
train_data = []
for train_graph in J['train_graph']:
    train_data.append([train_graph['before'], train_graph['pert'], train_graph['state'], train_graph['after']])

# all initial pts
initial_num = 2 ** n_node
initial_bin = []
for i in range(initial_num):
    initial_bin.append(bin(i)[2:])
initial_state = []
for i in initial_bin:
    default = '0' * (n_node - len(i)) + i
    initial_state.append([int(s) for s in default])

def binary(num, num2):
    initial_bin = []
    for i in range(num):
        initial_bin.append(bin(i)[2:])
    initial_state = []
    for i in initial_bin:
        default = '0' * (num2 - len(i)) + i
        initial_state.append([int(s) for s in default])
    return (initial_state)

if n_node > 10:
    cpu_num = 10
cpu_num = 2  ## have to smaller than node size, json simulation condition!
quo, re = divmod(n_node, cpu_num)
repeat = binary(2 ** cpu_num, cpu_num)
if re != 0:
    last = binary(2 ** re, re)
initial_list = []
for i in range(sample_size):
    initial_state = []
    for i in range(quo):
        r = repeat[np.random.choice(2 ** cpu_num, 1)[0]]
        initial_state.extend(r)
    if re != 0:
        L = last[np.random.choice(2 ** re, 1)[0]]
        initial_state.extend(L)
    initial_list.append(initial_state)


def get_attractor(Weight, initial_state, state_traj, pert_node=[], state=[]):
    # transition from t to t+1
    x = initial_state[:]
    if pert_node:
        for i, j in zip(pert_node, state):
            x[int(i)] = int(j)
    x.append(1)  # add basal level position
    x_array = np.array(x)
    traj_array = np.dot(Weight, x_array.transpose())

    traj_array[traj_array > 0] = int(1)
    for i in range(len(traj_array)):
        if traj_array[i] == 0:
            traj_array[i] = x_array[i]
    traj_array[traj_array < 0] = int(0)
    traj_list = [int(x) for x in traj_array]
    if pert_node:
        for i, j in zip(pert_node, state):
            traj_list[int(i)] = int(j)
    state_traj.append(initial_state)
    if traj_list == initial_state:
        return traj_list
    elif traj_list in state_traj:
        return state_traj[state_traj.index(traj_list):]
        # return traj_list
    else:
        return get_attractor(Weight=Weight, initial_state=traj_list, state_traj=state_traj, pert_node=pert_node,
                                  state=state)


### Model
def _func(w, initial_attractor, pert=[], state=[]):
    # make weight matrix
    W_matrix = np.zeros((n_node, n_node))
    for (enum, e) in enumerate(edge_list):
        source = node_list.index(e[0])
        target = node_list.index(e[2])
        W_matrix[target][source] = w[enum]
    basal_array = np.array(w[-n_node:])
    basal_array.shape = (n_node, 1)
    W_matrix = np.hstack((W_matrix, basal_array))
    state_traj = []
    result = get_attractor(Weight=W_matrix, initial_state=list(initial_attractor), state_traj=state_traj,
                                pert_node=pert, state=state)
    return result

f=open('output_all_noise3.txt','w')
# Reimplement the virtual method that defines the objective function
w = [-1,-2,-2,2,-2,2,1,1,-4,-4]
for i0 in range(-2,1):
    w[0]=i0
    for i1 in range(-2, 1):
        w[1] = i1
        for i2 in range(0, 3):
            w[2] = i2
            for i3 in range(-2, 1):
                w[3] = i3
                for i4 in range(0, 3):
                    w[4] = i4
                    for i5 in range(0, 3):
                        w[5] = i5
                        for i6 in range(-2, 3):
                            w[6] = i6
                            for i7 in range(-2, 3):
                                w[7] = i7
                                for i8 in range(-2, 3):
                                    w[8] = i8
                                    for i9 in range(-2,3):
                                        w[9] = i9
                                # get unique attractor
                                        attractor_set = []
                                        attractor_dic = {}
                                        for initial in initial_list:
                                            result_att = _func(w, initial, pert=list(input_dic.keys()), state=list(input_dic.values()))
                                            if type(result_att[0]) is int:
                                                if tuple(result_att) in attractor_dic:
                                                    continue
                                                attractor_dic[tuple(result_att)] = 0
                                                attractor_set.append(result_att)
                                            elif type(result_att[0]) is list:
                                                if tuple(result_att[0]) in attractor_dic:
                                                    continue
                                                for n, att in enumerate(result_att):
                                                    attractor_dic[tuple(att)] = 0
                                                    attractor_set.append(att)
                                # make reachability graph

                                        reachability = []
                                        for initial_att in attractor_set:
                                            for i in range(n_node):
                                                final_att = _func(w, initial_att, pert=list(input_dic.keys()) + [int(i)],
                                                                       state=list(input_dic.values()) + [0])
                                                if type(final_att[0]) is int:
                                                    reachability.append([initial_att, i, 0, final_att])
                                                elif type(final_att[0]) is list:
                                                    for att in final_att:
                                                        reachability.append([initial_att, i, 0, att])
                                                final_att = _func(w, initial_att, pert=list(input_dic.keys()) + [int(i)],
                                                                       state=list(input_dic.values()) + [1])
                                                if type(final_att[0]) is int:
                                                    reachability.append([initial_att, i, 1, final_att])
                                                elif type(final_att[0]) is list:
                                                    for att in final_att:
                                                        reachability.append([initial_att, i, 1, att])
                                        total_s = 0.0
                                        for state in train_data:
                                            obj = [[int(x) for x in list(state[0])], [int(x) for x in list(state[3])]]
                                            s = []
                                            for sim in reachability:
                                                if (int(sim[1]) == int(state[1])) & (int(sim[2]) == int(state[2])):
                                                    s1 = [(int(i) - int(j)) ** 2 for i, j in zip(obj[0], sim[0])]
                                                    s2 = [(int(i) - int(j)) ** 2 for i, j in zip(obj[1], sim[3])]
                                                    s.append(sum(s1) + sum(s2))
                                                    # print('score',sum(s1)+sum(s2), s)
                                                    # print('obj',obj)
                                                    # print('state',state)
                                                    # print('sim',sim)
                                            # print(min(s))
                                            total_s += int(min(s))
                                        # total_s = total_s + 0 * sum([abs(x) for x in w])
                                        total_s = total_s - len([x for x in w if x ==0])
                                        print(w, total_s)
                                        f.write('\t'.join([str(x) for x in w]))
                                        f.write('\t%s\n'%total_s)
                                # 2개와 score를 plot해보자