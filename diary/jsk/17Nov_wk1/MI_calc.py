### https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations

###Generating mutual information of every pair of antibodies

#function for calculating mutual information
def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi

#reading and processing raw data
raw = pd.read_csv('raw_data.csv',index_col=0)
c = [raw.columns.values[0], 'c']#data was separated by columns; the input conditions
q = [raw.columns.values[1:5], 'q']
r = [raw.columns.values[5:9], 'r']
i = [raw.columns.values[9:14], 'i']
a = [raw.columns.values[14:18], 'a']
s = [raw.columns.values[18:22], 's']

mraw=[[raw, '']]#the original data fitted to the relevant data structure
condcomb=[]
for it in range(2,5):#2, 3, and 4 combinations were made from the data
    condcomb = condcomb + [list(comb) for comb in combinations([q,r,i,a,s], it)]

for tempcomb in condcomb:#separate data tables were made for each raw dataset
    tempraw = [raw[np.concatenate(np.array(tempcomb).T[0].tolist())], ''.join(np.array(tempcomb).T[1].tolist())]
    mraw.append(tempraw)

#calculate pairwise mutual information for each dataset
for rawslice in mraw:
    rawdata = rawslice[0]
    inccond = rawslice[1]
    row_count = rawdata.shape[1]
    col_count = rawdata.shape[0]
    
    #generating pair combination of every pair of antibodies
    L = [k for k in range(col_count)]
    combination =[list(comb) for comb in combinations(L, 2)]
    
    #calculating mutual information and exporting them to .csv file
    result = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(rawdata.index),index=list(rawdata.index))
    for pair in combination:
        n = pair[0]
        m = pair[1]
        x = rawdata.iloc[[n]].values.flatten().tolist()
        y = rawdata.iloc[[m]].values.flatten().tolist()
        MI = calc_MI(x,y,row_count)
        result.iloc[n,m] = MI
    result.to_csv('MI_raw' + inccond + '.csv')






