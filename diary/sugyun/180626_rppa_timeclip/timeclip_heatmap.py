# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 14:34:36 2018

@author: ansu
"""

import pandas as pd

s1_b = pd.read_csv('s_0.1_biocarta.csv', index_col=0,encoding='cp1252')
s1_b.columns = ['path','alpha',1,6,24,96]
s2_b = pd.read_csv('s_0.2_biocarta.csv', index_col=0,encoding='cp1252')
s2_b.columns = ['path','alpha',1,6,24,96]
q1_b = pd.read_csv('q_0.1_biocarta.csv', index_col=0,encoding='cp1252')
q1_b.columns = ['path','alpha',1,6,24,96]
q2_b = pd.read_csv('q_0.2_biocarta.csv', index_col=0,encoding='cp1252')
q2_b.columns = ['path','alpha',1,6,24,96]
#p1_b = pd.read_csv('p_0.1_biocarta.csv', index_col=0,encoding='cp1252')
#p1_b.columns = ['path','alpha',0.5,1,3,6,24]
#p2_b = pd.read_csv('p_0.2_biocarta.csv', index_col=0,encoding='cp1252')
#p2_b.columns = ['path','alpha',0.5,1,3,6,24]

## alpha가 0.05 이하인 시그널 pathway 색출
## 결과를 보면 biocarta에서 유의미한 시그널들이 공통적으로 뽑히는 것을 관찰할 수 있다.
def df_filter(df):
    df = df.dropna()
    df = df[df['alpha']<0.05]
    return(df)


s1_b = df_filter(s1_b)
s2_b = df_filter(s2_b)
q1_b = df_filter(q1_b)
q2_b = df_filter(q2_b)
#p1_b = df_filter(p1_b)
#p2_b = df_filter(p2_b)


## heatmap 그리기
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()

s2_b_heat = s2_b.iloc[:,2:]
s2_b_heat.index = s2_b['path'].tolist()
q2_b_heat = q2_b.iloc[:,2:]
q2_b_heat.index = q2_b['path'].tolist()
ax_s = sns.clustermap(s2_b_heat,col_cluster=False,cmap="mako")
ax_q = sns.clustermap(q2_b_heat,col_cluster=False,cmap="mako")


