### https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations
def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi


fc = pd.read_csv('raw_data.csv',index_col=0)
q = fc[[1,2,3,4]]
# print(q)
# print(list(q.index))
row_count = q.shape[1]
col_count = q.shape[0]

L = [k for k in range(col_count)]
combination =[list(comb) for comb in combinations(L, 2)]

result = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(q.index),index=list(q.index))
a=10
for i in range(len(combination)):
    n = combination[i][0]
    m = combination[i][1]
    x = q.iloc[[n]].values.flatten().tolist()
    y = q.iloc[[m]].values.flatten().tolist()
    MI = calc_MI(x,y,row_count)
    # print(n,m,MI)
    result.ix[n,m] = MI
# print(result)
result.to_csv("MI_q.csv")
######################################################################
fc = pd.read_csv('raw_data.csv',index_col=0)
s = fc[[18,19,20,21]]
# print(s)
# print(list(s.index))
row_count = s.shape[1]
col_count = s.shape[0]

L = [k for k in range(col_count)]
combination =[list(comb) for comb in combinations(L, 2)]

result = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(s.index),index=list(s.index))
a=10
for i in range(len(combination)):
    n = combination[i][0]
    m = combination[i][1]
    x = s.iloc[[n]].values.flatten().tolist()
    y = s.iloc[[m]].values.flatten().tolist()
    MI = calc_MI(x,y,row_count)
    # print(n,m,MI)
    result.ix[n,m] = MI
# print(result)
result.to_csv("MI_s.csv")