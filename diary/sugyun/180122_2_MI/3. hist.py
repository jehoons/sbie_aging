import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt


fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flatten()

fc = pd.read_csv('MI_s.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
x=x.values.flatten().tolist()
ax0.hist(x, bins=50,range=(min(x),max(x)))
ax0.set_title('senescence')

fc = pd.read_csv('MI_q.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
# print(x)
x = x.values.flatten().tolist()
ax1.hist(x, bins=50,range=(min(x),max(x)))
ax1.set_title('quiescence')

fc = pd.read_csv('MI_apop.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
x=x.values.flatten().tolist()
ax2.hist(x, bins=50,range=(min(x),max(x)))
ax2.set_title('apop')

fc = pd.read_csv('MI_p.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
x=x.values.flatten().tolist()
ax3.hist(x, bins=50,range=(min(x),max(x)))
ax3.set_title('proliferation')

plt.show()