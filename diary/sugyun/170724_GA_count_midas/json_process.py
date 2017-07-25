import json
import os
dic = {}
dic['simulation_condition'] = {'sample_size':1000,
                               'repeat':500,
                               'lam':0.1,
                               'step_limit':100}
# dic['input_condition'] = [{'node':'x1', 'state':1}]
dic['input_condition'] = []

curent_path=os.curdir
os.chdir(curent_path+'/data')

dic['PKN'] = []
f= open('default_edge.csv','r')
for n, line in enumerate(f):
    if n == 0:
        continue
    else:
        line_list = line.strip().split((' '))
        dic['PKN'].append({'source': line_list[0], 'interaction': line_list[1], 'target': line_list[2]})

dic["Node"] = []
f= open('default_node.csv','r')
for n, line in enumerate(f):
    if n == 0:
        continue
    else:
        line_list = line.strip().split((' '))
        dic['Node'].append({'name':line_list[0]})

n_node = len(dic['Node'])
dic["train_graph"] = []
f= open('train_graph_midas.csv','r')
_dic ={}
for n, line in enumerate(f):
    if n == 0:
        continue
    else:
        line_list = line.strip().split((','))
        if line_list[0] in _dic:
            _dic[line_list[0]] =_dic[line_list[0]] + line_list[1:]
        else:
            _dic[line_list[0]]=line_list[1:]

for i in _dic:
    stack=_dic[i]
    pert_index = [i for i,j in enumerate(stack[:n_node]) if j !='0']
    pert_state = []
    for i in pert_index:
        pert_state.append(stack[i])
    stack = stack[n_node:]
    if ''.join(stack[:n_node]) == '0'*n_node:
        before = ''.join(stack[n_node:2*n_node])
        after = ''.join(stack[-n_node:])
    elif ''.join(stack[:n_node]) == '1'*n_node:
        after = ''.join(stack[n_node:2 * n_node])
        before = ''.join(stack[-n_node:])
    dic["train_graph"].append({'before':before,'pert': pert_index[0],'state':int(pert_state[0]),'after':after})

    # for j in _dic["train_graph"]
# dic["train_graph"] = [{'before':'0010','pert': 0,'state':0,'after':'0110'},
#                       {'before':'1010','pert': 0,'state':0,'after':'0110'},
#                       {'before':'0010','pert': 1,'state':0,'after':'0010'},
#                       {'before':'1100','pert': 1,'state':0,'after':'1000'},
#                       {'before':'1011','pert': 1,'state':0,'after':'1011'},
#                       {'before':'1011','pert': 1,'state':0,'after':'0010'},
#                       {'before':'1011','pert': 1,'state':0,'after':'1000'},
#                       {'before':'1010','pert': 1,'state':0,'after':'1010'},
#                       {'before':'0010','pert': 2,'state':0,'after':'0101'},
#                       {'before':'1100','pert': 2,'state':0,'after':'0101'},
#                       {'before':'1011','pert': 2,'state':0,'after':'0101'},
#                       {'before':'1010','pert': 2,'state':0,'after':'0101'},
#                       {'before':'0010','pert': 3,'state':0,'after':'1010'},
#                       {'before':'1010','pert': 3,'state':0,'after':'1010'}]
os.chdir('..')
with open('data.json', 'w') as outfile:
    # json.dump(data, outfile)
    json.dump(dic, outfile)