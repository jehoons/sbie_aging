import numpy as np
import pandas as pd
# import phenotype_detector as phd
# import attractor_phenotype as ap

def phenotype_detector(file, attractor_deter_list, logic_table):
    dic_result = attractor(file, attractor_deter_list)
    for i in logic_table:
        index = []
        check_series = ''
        for n, m in enumerate(i):
            if m != 'none':
                index.append(n)
                if len(m) < 2:
                    check_series = check_series + m
        index = index[:-1]
        for j in dic_result:
            ch = dic_result[j][4]
            target_series = ''
            for n in index:
                target_series = target_series + ch[n]
            if target_series == check_series:
                dic_result[j].append(i[-1])

    g = open(outputfile1, 'w')
    g.write('LowNuEx    IGF1    type    ratio   attractor   phenotyppe\n')
    for i in dic_result:
        for j in dic_result[i]:
            g.write("%s\t" %(j))
        g.write('\n')
    print('file is saved!')
    return 0


def attractor(file, attractor_deter_list):
    dic_result = {}
    for n, line in enumerate(file):
        line = line.strip()
        if not line:
            continue
        elif line.startswith('input_condition'):
            lownuex = line.split()[2][0:-1]  # LowNuEx
            IGF1 = line.split()[4][0:-1]  # IGF1
            continue
        elif line.startswith('type'):
            result_line = []
            result_line.append(lownuex)  # LowNuEx
            result_line.append(IGF1)  # IGF1
            result_line.append(line.split()[1])  # type
            result_line.append(line.split()[3])  # ratio
            continue
        elif line.startswith('sign'):
            protein_name_list = line.split(',')
            protein_name_list = protein_name_list[1:]
            index = []
            for i in attractor_deter_list:
                index.append(protein_name_list.index(i))  # get the index of attractor_deter_protein
            continue
        list_line = line.split(',')
        list_line = list_line[1:]
        if ('0' in list_line or '1' in list_line):
            state_rev = ''
            for i in index:
                state_rev += list_line[i]
            dic_result[n] = result_line +[state_rev]
    return dic_result


def test_c1(): 
    # Original code is written by 안수균 
    # Original Readme: 
    # Open the "sim_result_analysis-v2.py" and run.
    # This python file imports a module which is "phenotype_detector.py".
    # In addition, "phenotype_detector.py" imports a module which is "attractor_phenotype.py".
    # you have to put the file, gene symbol list, and the logic to get the phenotype of the attractor.
    # The result will be provided as txt file.    
    # by A.S.K

    f = open(inputfile1,'r')
    attractor_deter_list = ['caspase3','RB1','E2F','LC','mTOR']  # RB1 => pRB1 in log_v1.txt
    logic_table = np.array([
                            ['0', '0', '1', 'none', '1', 'Pro'],
                            ['1', 'none', 'none', 'none', 'none', 'Apop'],
                            ['0', '1', '0', 'none', '0', 'Qui'],
                            ['0', '1', '0', '1', '1', 'Qui'],
                            ['0', '1', '0', '0', '1', 'Sen']
                            ])

    result = phenotype_detector(f,attractor_deter_list,logic_table)  


inputfile1 = 'c/b5-simul-log_v2.txt'

outputfile1 = 'c/b5-updated1-boolsim-results-summary_phenotype_v2.txt'
