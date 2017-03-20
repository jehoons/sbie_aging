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

# output 
file_b1 = dirname(table_s1.__file__) + '/b/b1-%s-boolsim-inputs.json'
file_b2 = dirname(table_s1.__file__) + '/b/b2-%s-boolsim-results.json'
file_b3 = dirname(table_s1.__file__) + '/b/b3-%s-boolsim-results-summary.csv'


def test_simul():

    model_list = [(file_a3, 'initial'), (file_a4, 'updated1')]
    for file, keyword in model_list: 
        simulate(infile=file, keyword=keyword)
        result_summary(keyword=keyword)


def simulate(infile=None, keyword=None):
    # if exists(file_b1):
    #     return
    with open(infile, 'r') as fin: 
        lines = fin.readlines()

    nodes = [] 
    for l0 in lines: 
        words = l0.split('*=')
        nodes.append(words[0].strip())

    modeltext = ''
    for node in nodes: 
        modeltext+= '%s = Random\n' % node

    allnodes = "".join(lines).replace('and','').replace('or',''
        ).replace('not','').replace('(',' ').replace(')',' '
        ).replace('*=',' ').replace('\n','').split()

    inputnodes = [ x for x in set(allnodes) - set(nodes)]

    modeltext += "".join(lines)

    # make input configurations

    combi = [ {inputnodes[0]:i, inputnodes[1]:j}
        for i,j in itertools.product( (False,True), 
            repeat=len(inputnodes) )]
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
            'samples': 1000, 
            'steps': 500, 
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


def result_summary(keyword=None):
    datafile = file_b2%keyword
    print('inputfile:', datafile)

    with open(datafile, 'r') as fobj_c:
        data = json.load(fobj_c)

    df0 = pd.DataFrame([], columns=['IGF1','LowNuEx','attr_type','ratio','value'])

    i = 0 
    for d in data: 
        input_condition = d['input_condition']
        attrs = d['attractors']
        print('input_condition:', input_condition)
        for attr in attrs: 
            df0.loc[i, 'attr_type'] = attrs[attr]['type']
            df0.loc[i, 'ratio'] = attrs[attr]['ratio']
            df0.loc[i, 'value'] = attr
            df0.loc[i, 'IGF1'] = input_condition['IGF1']
            df0.loc[i, 'LowNuEx'] = input_condition['LowNuEx']

            print('type',attrs[attr]['type'],'ratio',attrs[attr]['ratio'])
            # print('value', attrs[attr]['value'])

            value = attrs[attr]['value']
            attr_type = attrs[attr]['type']

            if attr_type == 'point': 
                labels = d['labels']
                binstr = d['state_key'][value]
                print(",".join(labels))
                print (",".join([b for b in binstr]))
            elif attr_type == 'cyclic': 
                print(",".join(labels))
                for c0 in value: 
                    binstr = d['state_key'][c0]
                    print (",".join([b for b in binstr]))

            print()
            i+=1

            # set_trace() ; break 
    summary_output = file_b3%keyword
    df0.to_csv(summary_output)
