import numpy as np
import pandas as pd
from statistics import variance

#data import
sene_data = pd.read_csv('s_sym_normlog.csv', index_col=0)
qui_data = pd.read_csv('q_sym_normlog.csv', index_col=0)
pro_data = pd.read_csv('p_sym_normlog.csv', index_col=0)
pct=.2#percentile; change this

for data, cond in zip([sene_data, qui_data, pro_data], ['s', 'q', 'p']):
    datalen = len(data)
    dataidxs = list(data.index)
    var_data = list(zip(data.apply(variance, axis=1), np.arange(datalen)))
    sorted_var_data = sorted(var_data, reverse=True)
    indices = [idx for _, idx in sorted_var_data[:int(datalen*pct)]]
    tmpdata = np.array(data)[indices]
    pd.DataFrame(tmpdata, index=np.array(dataidxs)[indices], columns=data.columns).to_csv(cond+'_'+str(pct)+'_var_cutoff.csv')














