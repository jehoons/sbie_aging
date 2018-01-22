### https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations

diff = pd.read_csv('MI_diff.csv',index_col=0)
col_count = diff.shape[1]
diff = diff.abs()
diff_complete = diff.add(pd.DataFrame.transpose(diff), fill_value=0)
diff_complete.to_csv("MI_diff_abs_total.csv")
diff_complete = diff_complete[diff_complete>0.3]
diff_complete['sum'] = diff_complete.sum(axis=0)
diff_complete=diff_complete.sort_values(by=['sum'],ascending=False)
diff_complete.to_csv("MI_diff_target.csv")

# target = diff_complete.sum(axis=1)
# target.to_csv("MI_diff_target.csv")
# target =
#
#
#
# diff_abs = abs(s-q)
# diff2 = diff[diff_abs>0.3]
# diff2.to_csv("MI_diff.csv")