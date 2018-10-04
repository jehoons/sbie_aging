# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 02:01:13 2018

@author: ansu
"""

import pandas as pd
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()

data = pd.read_csv('auc_heatmap.csv', index_col = 0)
data = data[['Quiescence','Senescence','Proliferation']]

ax = sns.heatmap(data, center=1)

from matplotlib.colors import LinearSegmentedColormap
data = pd.read_csv('auc_heatmap_discrete.csv', index_col = 0)
data = data[['Quiescence','Senescence','Proliferation']]
colors = ["black", "lightgray"] 
cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))

ax = sns.heatmap(data, cmap=cmap,cbar_kws={"shrink": .5})
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.25,0.75])
colorbar.set_ticklabels(['0', '1'])