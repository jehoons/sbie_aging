"""
아모래의 rppa fold change 데이터를 이용하여 
senescence와 quiescence 각 컨디션에 대하여
모든 단백질 조합의 Mutual information을 구함

""" 


import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations

# mi 계산 함수
def calc_MI(x, y, bins):
    c_xy = np.histogram2d(x, y, bins)[0]
    mi = mutual_info_score(None, None, contingency=c_xy)
    return mi
# 결과 함수
def result_fnc(combination, fc, matrix):
    result = matrix
    for i in range(len(combination)):
        print(i/len(combination))
        n = combination[i][0]
        m = combination[i][1]
        x = fc.iloc[[n]].values.flatten().tolist()
        y = fc.iloc[[m]].values.flatten().tolist()
        MI = calc_MI(x,y,row_count)
        result.iloc[n,m] = MI
    return(result)
    

# 파일 입력 및 변수 생성
fc = pd.read_csv('1. fold_change.csv',index_col=0)
q = fc.iloc[:,[1,2,3,4]]
s = fc.iloc[:,[18,19,20,21]]
p = fc.iloc[:,[9,10,11,12,13]]
a = fc.iloc[:,[14,15,16,17]]
row_count = q.shape[1]
col_count = q.shape[0]
combination =[list(comb) for comb in combinations([k for k in range(col_count)], 2)]
matrix = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(q.index),index=list(q.index))
# 결과
result = result_fnc(combination, q, matrix)
result.to_csv("1. MI_q.csv")

row_count = s.shape[1]
col_count = s.shape[0]
combination =[list(comb) for comb in combinations([k for k in range(col_count)], 2)]
matrix = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(s.index),index=list(s.index))
#print(result)
result = result_fnc(combination, s, matrix)
result.to_csv("1. MI_s.csv")

row_count = p.shape[1]
col_count = p.shape[0]
combination =[list(comb) for comb in combinations([k for k in range(col_count)], 2)]
matrix = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(p.index),index=list(p.index))
#print(result)
result = result_fnc(combination, p, matrix)
result.to_csv("1. MI_p.csv")

row_count = a.shape[1]
col_count = a.shape[0]
combination =[list(comb) for comb in combinations([k for k in range(col_count)], 2)]
matrix = pd.DataFrame(np.zeros((col_count,col_count)),columns=list(a.index),index=list(a.index))
#print(result)
result = result_fnc(combination, a, matrix)
result.to_csv("1. MI_a.csv")
#print(result)