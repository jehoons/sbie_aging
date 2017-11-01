from positive_pruning_v3 import posprune
#from negative_pruning_v3 import negprune
import pandas as pd
from itertools import combinations

#importing relevant network-related data
agingnetwork = pd.read_csv('Aging_network.csv')
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)

condcomb=[['']]#include the original MI_raw

for it in range(2,5):#2, 3, and 4 combinations were made from previous data preparation
    condcomb = condcomb + [list(comb) for comb in combinations(['q','r','i','a','s'], it)]

mirawlist = []
for tempcomb in condcomb:#import each mutual information dataset
    mirawname = ''.join(tempcomb)
    tempmiraw = pd.read_csv('MI_raw' + mirawname + '.csv', index_col=0)
    mirawlist.append([tempmiraw, mirawname])#bind the imported MI_raw file with its indicative mirawname

for tempmiraw, mirawname in mirawlist:#perform positive and negative pruning and export the results
    agnetp, cutoffp = posprune(tempmiraw, agingnetwork, abmap, .5)
    #agnetn, cutoffn = negprune(tempmiraw, agingnetwork, abmap, .5)
    agnetp.to_csv('agingnetwork' + mirawname + str(int(cutoffp*10)) + 'p.sif', sep='\t', index=None, header=None)
    #agnetn.to_csv('agingnetwork' + mirawname + str(int(cutoffn*10)) + 'n.sif', sep='\t', index=None, header=None)




