## attractor fnc
import numpy as np
def get_attractor(Weight, initial_input, state_traj, pert_node = -1):
    # transition from t to t+1
    x = initial_input[:]
    if pert_node != -1:
        x[pert_node] = 0
    # x[1] = 0
    # x[2] = 0
    # x[3] = 0
    x.append(1)# add basal level position
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
        return get_attractor(Weight=Weight,initial_input=traj_list, state_traj=state_traj, pert_node =pert_node)

def basin_size(attractor_dic):
    uniq_attractor = {}
    cyc_att_set = set()
    for i in attractor_dic.values():
        # point attractor
        if type(i[0]) is int:
            i_str = ''.join([str(x) for x in i])
            if i_str in uniq_attractor:
                uniq_attractor[i_str] += 1
            else:
                uniq_attractor[i_str] = 1
            continue
        # cyclic attractor
        elif type(i[0]) is list:
            # 이름을 만듬
            att_name = ''
            for n_att, k in enumerate(i):
                k_str = ''.join([str(x) for x in k])
                # 하나라도 state가 같으면 같은 cyclic
                if k_str in cyc_att_set:
                    for l in uniq_attractor:
                        if k_str in l:
                            if n_att == 0:
                                uniq_attractor[l] += 1
                    continue
                # 없으면 새로운 cyclic att 이름 만들기
                elif k_str not in cyc_att_set:
                    cyc_att_set.add(k_str)
                    att_name = att_name + k_str
            if att_name == '':
                continue
            if att_name not in uniq_attractor:
                uniq_attractor[att_name] = 1
    return uniq_attractor