import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations
import math

ratio = pd.read_csv('MI_ratio.csv',index_col=0)
thres = 4
ratio1 = ratio[ratio>=thres]
ratio1_rev = ratio1.add(pd.DataFrame.transpose(ratio1), fill_value=0)
ratio1_rev = ratio1_rev.apply(np.log2)

ratio1_rev.to_csv("MI_ratio_up.csv")
ratio2 = ratio[ratio<=1/thres]
ratio2_rev = ratio2.add(pd.DataFrame.transpose(ratio2), fill_value=0)
ratio2_rev = ratio2_rev.apply(np.log2)
ratio2_rev.to_csv("MI_ratio_down.csv")



new = ratio1_rev.add(ratio2_rev, fill_value=0)
# new.to_csv("MI_ratio_total.csv")
new = new.apply(abs)
new['sum'] = new.sum(axis=0)
new.to_csv("MI_ratio_total_abs.csv")

