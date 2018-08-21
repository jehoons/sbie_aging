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

#redefining variables
nodelist = list(conddata.columns)


#manual stuff; just temporary
optfile = open('output_fittedparams.txt', 'r')
idc = 0
for line in optfile:
    if idc == 1:
        w = [float(i) for i in line.strip().split('\t')]
    if idc == 2:
        b = [float(j) for j in line.strip().split('\t')]
    idc = 0
    if line == 'Mean Weights:\n':
        idc = 1
    if line == 'Mean Basal Activity:\n':
        idc = 2

ws=[w]
bs=[b]
'''

#again manual stuff; just temporary
optfile = open('output_fittedparams.txt', 'r')
idc = 0
ws = []
bs = []
for line in optfile:
    if idc == 1:
        ws.append([int(i) for i in line.strip().split('\t')])
    if idc == 2:
        bs.append([int(j) for j in line.strip().split('\t')])
    idc = 0
    if line == 'Weights:\n':
        idc = 1
    if line == 'Basal Activities:\n':
        idc = 2

#single indegree nodes; the single indegree node's weight and basal activity waters down to a single value; adjust accordingly
indeglist = {}
for node in nodelist:
    indeglist[node] = 0
for link in netdata:
    indeglist[link[2]] += 1
singleindeg = [indeg for indeg in indeglist if indeglist[indeg] == 1]

#the stuff im working on; about the adjustment stuff described above
indegidx = list(np.array(netdata).transpose()[2])
phenlist = []
for w, b in zip(ws, bs):
    phens = []
    for sing in singleindeg:
        sw = w[indegidx.index(sing)]
        if netdata[indegidx.index(sing)][1] == 'inhibit':
            sw = -sw
        sb = b[nodelist.index(sing)]
        phen = sw + sb
        phens.append(phen)
    phens = np.array(phens)
    phens[phens > 0] = 1
    phens[phens <= 0] = 0
    phenlist.append(list(phens))

'''



###
for w, b in zip(ws, bs):
    #utilize _func in OptWeightandBasal for attractor calculation
    attlist = []
    wab = OptWeightandBasal(dim = dim, i_dim = dim, netdata = netdata, conddata = conddata)
    for condno in range(wab.num_conds):#for all of the conditions presented calculate the cost function and add all
        condno += 1
        inistates = wab.conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
        attractors = wab._func(w, b, inistates)#calculated states form the Boolean model
        attlist.append(attractors)
    
    fitlist = []
    for atts in attlist:
        fitness = 0.0
        for condno in range(wab.num_conds):#for all of the conditions presented calculate the cost function and add all
            condno += 1
            inistates = conddata.loc['INPUT: ' + str(condno)].tolist()#initial states given
            objstates = np.array(conddata.loc['OUTPUT: ' + str(condno)].tolist())#objective states
            objboolmask = np.array([state != 9 for state in objstates])#Boolean mask for selecting only the states required for cost function calculation
            objphen = list(objstates[objboolmask])
            tmpfitlist = []#list of fitness value for all of the states within a cyclic attractor; only one of if a point attractor
            for att in atts:
                att = np.array(att)
                predphen = list(att[objboolmask])
                tmpfit = sum((objphen[idx] - predphen[idx])**2 for idx in range(len(objphen)))#fitness function calulated for this condition and added
                tmpfitlist.append(tmpfit)
            fitness += np.mean(np.array(tmpfitlist))#fitness is the mean of fitness of all the states in a cyclic attractor
        fitlist.append(fitness)
    
    print('w')
    print(w)
    print('b')
    print(b)
    print('attlist')
    print(attlist)
    print('fitlist')
    print(fitlist)






