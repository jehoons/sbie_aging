import json
from os.path import exists
import os
from sbie_aging.results import table_s1
from os.path import dirname,basename
from ipdb import set_trace
import pandas as pd 
import itertools 


# input 
file_a3 = dirname(table_s1.__file__) + '/a/a3-logical-model-simple-rule.txt'
file_a4 = dirname(table_s1.__file__) + '/a/a4-logical-model-updated-1.txt'
file_a5 = dirname(table_s1.__file__) + '/a/a4-logical-model-updated-2.txt'
file_a6 = dirname(table_s1.__file__) + '/a/a4-logical-model-updated-3.txt'
file_a7 = dirname(table_s1.__file__) + '/a/a4-logical-model-updated-4.txt'

# output 
file_b1 = dirname(table_s1.__file__) + '/b/b1-%s-boolsim-inputs.json'
file_b2 = dirname(table_s1.__file__) + '/b/b2-%s-boolsim-results.json'
file_b3 = dirname(table_s1.__file__) + '/b/b3-%s-boolsim-results-summary.csv'
file_b4 = dirname(table_s1.__file__) + '/b/b4-%s-state-simmat.csv'


def test_simul():

    model_list = [(file_a3, 'initial', ['IGF1','LowNuEx']),
                  (file_a4, 'updated1', ['IGF1','LowNuEx']),
                  (file_a5, 'updated2', ['IGF1','LowNuEx','DNAdamage']),
                  (file_a6, 'updated3', ['IGF1','LowNuEx','DNAdamage']),
                  (file_a7, 'updated4', ['IGF1','LowNuEx','DNAdamage'])]

    for file, keyword,inputcols in model_list:
        simulate(infile=file, keyword=keyword)
        result_summary(keyword=keyword, inputcols=inputcols)


def simulate(infile=None, keyword=None):
    # if exists(file_b1):
    #     return
    with open(infile, 'r') as fin: 
        lines = fin.readlines()

    nodes = []
    # lines2 = []
    for l0 in lines:
        # if l0.startswith('#'):
        #     print('skipped: '+l0)
        #     continue
        words = l0.split('*=')
        nodes.append(words[0].strip())
        # lines2.append(lines)

    modeltext = ''
    for node in nodes: 
        modeltext+= '%s = Random\n' % node

    allnodes = "".join(lines).replace('and','').replace('or',''
        ).replace('not','').replace('(',' ').replace(')',' '
        ).replace('*=',' ').replace('\n','').split()

    inputnodes = [ x for x in set(allnodes) - set(nodes)]

    modeltext += "".join(lines)

    # make input configurations
    # set_trace()
    if len(inputnodes) == 3:
        combi = [ {inputnodes[0]:i, inputnodes[1]:j, inputnodes[2]:k}
                  for i, j, k in itertools.product((False, True),
                                                   repeat=len(inputnodes))]
    elif len(inputnodes) == 2:
        combi = [{inputnodes[0]: i, inputnodes[1]: j}
            for i, j in itertools.product((False, True),
                repeat=len(inputnodes))]
    else:
        print('unknown number of inputnodes')
        set_trace()
        assert False

    input_list = []
    for combi0 in combi:     
        # print(combi0)
        reverse_dict = {False:[], True:[]} 
        for k in combi0:
            v = combi0[k]
            reverse_dict[v].append(k)

        # print(reverse_dict)    
        input_list.append({
            'model': modeltext.split('\n'),
            'on_states': reverse_dict[True],
            'off_states': reverse_dict[False],
            'input_condition': combi0,
            'samples': 10000, 
            'steps': 150, 
            'debug': False
            })

    simul_infile = file_b1%keyword
    simul_outfile = file_b2%keyword

    with open(simul_infile, 'w') as foutstring: 
        json.dump(input_list, foutstring, indent=4)
    
    if not exists(simul_outfile): 
        os.system('python myengine.py %s %s' % (
            simul_infile, simul_outfile) )
    else: 
        # print('already exists: %s' % simul_outfile)
        pass


def result_summary(keyword=None, inputcols=None):
    datafile = file_b2 % keyword
    print('#inputfile:', datafile)

    with open(datafile, 'r') as fobj_c:
        data = json.load(fobj_c)

    # df0 = pd.DataFrame([], columns=['IGF1','LowNuEx','DNA_damage','attr_type','ratio','value','size'])
    df0 = pd.DataFrame([], columns=inputcols+['attr_type', 'ratio', 'value', 'size'])

    i = 0 
    for d in data: 
        input_condition = d['input_condition']
        attrs = d['attractors']
        print('#input_condition:', input_condition)
        for attr in attrs: 
            df0.loc[i, 'attr_type'] = attrs[attr]['type']
            df0.loc[i, 'ratio'] = attrs[attr]['ratio']
            df0.loc[i, 'value'] = attr

            for inputcols0 in inputcols:
                df0.loc[i, inputcols0] = input_condition[inputcols0]

            print('#type',attrs[attr]['type'],'ratio',attrs[attr]['ratio'])
            # print('value', attrs[attr]['value'])

            value = attrs[attr]['value']
            attr_type = attrs[attr]['type']
            labels = ['sign'] + d['labels']

            # set_trace()

            if attr_type == 'point':
                binstr = d['state_key'][value]
                print(",".join(labels))
                print (",".join([value]+[b for b in binstr]))
            elif attr_type == 'cyclic':
                print(",".join(labels))
                for cyckey in value: 
                    binstr = d['state_key'][cyckey]
                    print (",".join([cyckey]+[b for b in binstr]))

            if attrs[attr]['type'] == 'point':
                df0.loc[i, 'size'] = 1
            elif attrs[attr]['type'] == 'cyclic':
                df0.loc[i, 'size'] = len(attrs[attr]['value'])

            print()
            i+=1

    summary_output = file_b3%keyword
    df0.to_csv(summary_output)

    d['state_key']
    import numpy as np 
    mat = np.zeros([])
    state_key = d['state_key']
    mat = np.zeros([len(state_key), len(state_key)])
    for i,s1 in enumerate(state_key):
        for j,s2 in enumerate(state_key):
            vec1 = state_key[s1]
            vec2 = state_key[s2]
            x1 = np.array([int(v) for v in vec1])
            x2 = np.array([int(v) for v in vec2])
            dist = sum( (x1 - x2)**2 ) 
            # print (dist)
            mat[i,j] = dist

    simmat = pd.DataFrame(mat, columns=state_key.keys(), index=state_key.keys())
    simmat.to_csv(file_b4%keyword)

