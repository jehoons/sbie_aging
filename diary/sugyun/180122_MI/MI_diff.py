### https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations

q = pd.read_csv('MI_q.csv',index_col=0)
s = pd.read_csv('MI_s.csv',index_col=0)
col_count = q.shape[1]
diff = pd.DataFrame(pd.DataFrame(np.zeros((col_count,col_count)),columns=list(q.index),index=list(q.index)))
diff = s-q
diff_abs = abs(s-q)
diff2 = diff[diff_abs>0.3]
diff2.to_csv("MI_diff.csv")
# print(diff)

q = q[q>0]
s = s[s>0]
# print(s)
ratio = pd.DataFrame(pd.DataFrame(np.zeros((col_count,col_count)),columns=list(q.index),index=list(q.index)))
ratio = s.divide(q)
ratio.to_csv("MI_ratio.csv")
# print(ratio)