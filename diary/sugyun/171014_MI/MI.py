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
# print(list(fc.index))
row_count = fc.shape[1]
col_count = fc.shape[0]

L = [k for k in range(col_count)]
combination =[list(comb) for comb in combinations(L, 2)]

result = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(fc.index),index=list(fc.index))
a=10
for i in range(len(combination)):
    n = combination[i][0]
    m = combination[i][1]
    x = fc.iloc[[n]].values.flatten().tolist()
    y = fc.iloc[[m]].values.flatten().tolist()
    MI = calc_MI(x,y,row_count)
    print(n,m,MI)
    result.ix[n,m] = MI
# print(result)
result.to_csv("MI_raw.csv")