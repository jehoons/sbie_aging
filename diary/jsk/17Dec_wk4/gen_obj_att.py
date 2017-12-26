import numpy as np
import pandas as pd

"""
This code generates sample input data based on template objective attractor set
The template contains input and output state of Boolean network with 0 and 1 in decided positions
and 9 in undecided positions. This code generates random input conditions by replacing 9 with 0 and 1 randomly.
The exported .csv file is used in the main.py.
"""
sampleno = 100#number of sample inputs to generate per condition

#importing template data; 0 and 1 means state variable, 9 means undefined
data = pd.read_csv('obj_attractor_template.csv', index_col = 0)

gendatalist = []#list of generated sample input conditions
iptconds = [cond for cond in data.index if 'INPUT:' in cond]
for iptcond in iptconds:
    for sp in range(sampleno):
        tmpipt = data.loc[iptcond].tolist()
        for idx, nd in enumerate(tmpipt):
            if nd == 9:
                tmpipt[idx] = np.random.choice(2)
            else:
                pass
        optname = 'OUTPUT: ' + iptcond[7:]
        gendatalist.append([tmpipt, optname])

obj_attractor = pd.DataFrame(index = data.columns, dtype = int)#dataframe with sampled input with appropriate output
for idx, iptopt in enumerate(gendatalist):
    idx += 1
    obj_attractor['INPUT: ' + str(idx)] = pd.Series(iptopt[0], index = obj_attractor.index)
    obj_attractor['OUTPUT: ' + str(idx)] = data.loc[iptopt[1]]
obj_attractor = obj_attractor.transpose()

#export generated input and output conditions
obj_attractor.to_csv('obj_attractor.csv')
