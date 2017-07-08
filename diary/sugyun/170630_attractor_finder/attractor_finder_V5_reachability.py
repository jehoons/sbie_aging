import numpy as np
import numpy as np

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
print('initial condition: ', initial_state)

## attractor fnc
def get_attractor(Weight, initial_input, state_traj):
    # transition from t to t+1
    x = initial_input[:]
    # x[0] = 0
    # x[1] = 0
    # x[2] = 0
    x[3] = 0
    x.append(1)# add basal level position
    x_array = np.array(x)
    traj_array = np.dot(Weight, x_array.transpose())

    traj_array[traj_array > 0] = int(1)
    traj_array[traj_array <= 0] = int(0)
    traj_list = [int(x) for x in traj_array]
    # traj_list[0] = 0
    # traj_list[1] = 0
    # traj_list[2] = 0
    traj_list[3] = 0


    state_traj.append(initial_input)
    if traj_list == initial_input:
        return traj_list
    elif traj_list in state_traj:
        return state_traj[state_traj.index(traj_list):]
        # return traj_list
    else:
        return get_attractor(Weight=Weight,initial_input=traj_list, state_traj=state_traj)

## get attractor
attractor_dic = {}
for initial in initial_state:
    state_traj = []
    # W는 np array, initial와 state_traj는 리스트
    attractor = get_attractor(Weight=W,initial_input=initial,state_traj=state_traj)
    attractor_dic[tuple(initial)] = attractor
# f = open("attractor_pertx3.txt",'w')
# for i in attractor_dic:
#     init = ''.join([str(int(x)) for x in i])
#     att = ''.join([str(int(x)) for x in attractor_dic[i]])
#     f.write("%s\t%s\n"%(init,att))
for i in attractor_dic:
    print('attr',i,attractor_dic[i])

## basin 크기 구하기
uniq_attractor = {}
cyc_att_set = set()
cycle = False
for i in attractor_dic.values():
    # point attractor
    if type(i[0]) is int:
        i_str = ''.join([str(x) for x in i])
        if i_str in uniq_attractor:
            uniq_attractor[i_str] = uniq_attractor[i_str]+1
        else:
            uniq_attractor[i_str] = 1
        continue
    # cyclic attractor
    elif type(i[0]) is list:
        # 이름을 만듬
        att_name = ''
        for n_att,k in enumerate(i):
            k_str = ''.join([str(x) for x in k])
            # 하나라도 state가 같으면 같은 cyclic
            if k_str in cyc_att_set:
                for l in uniq_attractor:
                    if k_str in l:
                        if n_att == 0:
                            uniq_attractor[l] = uniq_attractor[l] + 1
                continue
            # 없으면 새로운 cyclic att 이름 만들기
            elif k_str not in cyc_att_set:
                cyc_att_set.add(k_str)
                att_name = att_name+k_str+'/'
        if att_name == '':
            continue
        if att_name not in uniq_attractor:
            uniq_attractor[att_name] = 1
print(uniq_attractor)