import numpy as np
import pandas as pd
from optweightandbasal import OptWeightandBasal

#import network and objective network conditions with inputs
conddata = pd.read_csv('obj_attractor.csv', index_col = 0)#pd.DataFrame with input and output condition information
netdata = []#list with network information
netfile = open('PKN_19.sif', 'r')
for link in netfile:
    netdata.append(link.strip().split('\t'))
dim = len(netdata) + len(conddata.columns)#genetic algorithm parameter dimension; number of edges(weights) + number of nodes(basal activities)

#manual stuff; just temporary
w = [3,1,3,3,2,3,2,3,2,4,1,3,4,1,3,2,3,2,1,1,1,2,2,4,2,1,4,3,3,3,3,4,1,1,3,2,1,1,2,2,1,2,3,1,3,2,3,2,3,4,1,2,4,3,1,1,3,2,3,2,1,2,1,2,3,4,1,3,2,1,2,3,2,2,1,4,4,1]
b = [0,-1,-2,-2,0,0,0,-2,-3,-1,-2,-3,-2,-2,-2,-3,0,-1,-1,-1,0,-2,0,-3,-1,0,-1,0,-3,0,-1,-1,-2,-2,0,-1,-3,-1,-3,-1,-1,0,-2,-1,-2,0,-1,-3,0,-2]

#utilize _func in OptWeightandBasal for attractor calculation
attlist = []
wab = OptWeightandBasal(dim = dim, i_dim = dim, netdata = netdata, conddata = conddata)
for condno in range(wab.num_conds):#for all of the conditions presented calculate the cost function and add all
    condno += 1
    inistates = wab.conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
    attractors = wab._func(w, b, inistates)#calculated states form the Boolean model
    attlist.append(attractors)
print(attlist)
