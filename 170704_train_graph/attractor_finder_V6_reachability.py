import numpy as np
import numpy as np
from attractor_detect import get_attractor,basin_size
## edge & node input Data
file_edge = open('edge_info.txt', 'r')
file_node = open('node_info.txt', 'r')
node_list = []
edge_list = []
for line in file_node:
    node_list.append(line.strip())

for line in file_edge:
    edge_list.append(line.strip().split(' '))
#########################################################################
## make weight matrix
n_node = len(node_list)
W = np.zeros((n_node, n_node))
# edge weight 4개
w_list = [1, 1, -2, 2]
w_list = [-2, -2, 2, -2, 2]
# w_list = [1]*len(edge_list)
if len(w_list) != len(edge_list):
    print("check the number of link weight again")
    raise(IndexError)
# node 수와 같이 3개
basal_list = [0, -2, 1]
basal_list = [1, 1, -1, 1]
# basal_list = [1]*len(node_list)
if len(node_list) != len(basal_list):
    print("check the number of link weight again")
    raise(IndexError)
for (enum, e) in enumerate(edge_list):
    source = node_list.index(e[0])
    target = node_list.index(e[2])
    W[target][source] = w_list[enum]
basal_array = np.array( basal_list)
basal_array.shape = (n_node,1)
W = np.hstack((W, basal_array))
##################################################################################
## all initial combination
# 2에 11승
initial_num = 2**n_node
initial_bin = []
for i in range(initial_num):
    initial_bin.append(bin(i)[2:])
initial_state = []
for i in initial_bin:
    default = '0'*(n_node-len(i))+i
    initial_state.append([int(s) for s in default])
# print('initial condition: ', initial_state)
####################################################################################
## get attractor
attractor_dic = {}
for initial in initial_state:
    state_traj = []
    # W는 np array, initial와 state_traj는 리스트,pert_node는 int입니다.
    attractor = get_attractor(Weight=W,initial_input=initial,state_traj=state_traj, pert_node = -1)
    attractor_dic[tuple(initial)] = attractor
# for i in attractor_dic:
    # print('attr',i,attractor_dic[i])
################################################################
## basin 크기 구하기
basin = basin_size(attractor_dic=attractor_dic)
# print(basin)
initial_train = set(basin.keys())
ini_train_list = []
for s in initial_train:
    while s:
        ini_train_list.append(s[0:n_node])
        s = s[n_node:]
# print(ini_train_list)
##################################################################
## perturbation effect
g = open('train_graph.txt','w')
for i in range(n_node):
    for b in ini_train_list:
        initial = [int(x) for x in list(b)]
        state_traj = []
        attractor = get_attractor(Weight=W, initial_input=initial, state_traj=state_traj, pert_node=i)
        if type(attractor[0]) is int:
            print(i, ''.join([str(x) for x in initial]),''.join([str(x) for x in attractor]))
            g.write('%s\t%s\t%s\n'%(i, ''.join([str(x) for x in initial]),''.join([str(x) for x in attractor])))
        if type(attractor[0]) is list:
            for v in attractor:
                g.write('%s\t%s\t%s\n' % (i, ''.join([str(x) for x in initial]), ''.join([str(x) for x in v])))