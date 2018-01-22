import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flatten()

fc = pd.read_csv('MI_s.csv',index_col=0)
s = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
s=s.values.flatten().tolist()
s = [x for x in s if str(x) != 'nan']

fc = pd.read_csv('MI_q.csv',index_col=0)
q = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
q=q.values.flatten().tolist()
q = [x for x in q if str(x) != 'nan']

fc = pd.read_csv('MI_apop.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
a=x.values.flatten().tolist()
a = [x for x in a if str(x) != 'nan']

fc = pd.read_csv('MI_p.csv',index_col=0)
x = fc.where(np.triu(np.ones(fc.shape)).astype(np.bool))
p=x.values.flatten().tolist()
p = [x for x in p if str(x) != 'nan']

ax0.scatter(s,q)
ax1.scatter(s,a)
ax2.scatter(s,p)
ax3.scatter(p,a)
cor = pearsonr(s,q)
print(cor)
cor = pearsonr(s,a)
print(cor)
cor = pearsonr(s,p)
print(cor)
cor = pearsonr(q,a)
print(cor)
cor = pearsonr(q,p)
print(cor)
cor = pearsonr(a,p)
print(cor)


plt.show()
