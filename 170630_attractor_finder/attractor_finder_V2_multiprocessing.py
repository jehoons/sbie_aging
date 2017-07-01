import time
import numpy as np
from multiprocessing import Process
from multiprocessing import Pool

## edge & node input Data
file_edge = open('edge_info.txt', 'r')
file_node = open('node_info.txt', 'r')
node_list = []
edge_list = []
for line in file_node:
    node_list.append(line.strip())

for line in file_edge:
    edge_list.append(line.strip().split(' '))

## make weight matrix
n_node = len(node_list)
W = np.zeros((n_node, n_node))
# edge weight 4개
# w_list = [-10, 10, 2, 10, -9, 2, 3, 2, 5, 6, 9, 10, 8]
w_list = [1, 1, -2, 2]
# w_list = [1]*len(edge_list)

if len(w_list) != len(edge_list):
    print("check the number of link weight again")
    raise(IndexError)
# node 수와 같이 3개
# basal_list = [ -5, -3, -3, -9, -3, 4, -3, -4, -9, 7, -1]
basal_list = [0, -2, 1]
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

def get_attractor(Weight, initial_input, state_traj):
    # transition from t to t+1
    x = initial_input[:]
    x.append(1)# add basal level position
    x_array = np.array(x)
    traj_array = np.dot(Weight, x_array.transpose())

    traj_array[traj_array > 0] = int(1)
    traj_array[traj_array <= 0] = int(0)
    traj_list = list(traj_array)
    state_traj.append(initial_input)
    if traj_list == initial_input:
        return traj_list
    elif traj_list in state_traj:
        return state_traj[state_traj.index(traj_list):]
    else:
        return get_attractor(Weight=Weight,initial_input=traj_list, state_traj=state_traj)

# get attractor
# attractor_dic = {}
def f(num):
    initial = initial_state[num]
    state_traj = []
    # W는 np array, initial와 state_traj는 리스트
    attractor = get_attractor(Weight=W,initial_input=initial,state_traj=state_traj)
    return initial, attractor



if __name__=='__main__':
    start = time.time()
    p = Pool(4)
    result = p.map(f,range(2**n_node))
    attractor_dic = {}
    for i in result:
        attractor_dic[tuple(i[0])] = i[1]
    for i in attractor_dic:
        print(i, attractor_dic[i])
    end = time.time() - start
    print(end)









