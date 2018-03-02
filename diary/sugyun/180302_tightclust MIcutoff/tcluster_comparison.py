import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from sklearn.manifold import TSNE
from matplotlib import cm
from sklearn.manifold import TSNE
from matplotlib import cm
import seaborn as sns
import mygene as mg

# load data
result = pd.read_csv('result_with_sym.csv',index_col=0)

list_sym = result['symbol'].tolist()
rdm = np.random.choice(list_sym,100)


# data 수정
extra = result.pop('tight5')
extra = result.pop('tight4')
result = result.sort_values(by=['kclust'],ascending=True)

##각 클러스터의 각 tight clust 저장
#Cresult = result[result['tight10']>0]
#for i in range(4):
#    kc = Cresult[Cresult['kclust']==i]
#    for j in range(1,11):
#        tc = kc[kc['tight10']==j]
#        if tc.shape[0] == 0:
#            continue
#        tc[['kclust','tight10','entrez','symbol']].to_csv('C%s_T%s.csv'%(i,j))
genedf = result[result['tight10']>0]
gc1 = genedf[genedf['kclust']==1]
gc1[['tight10','kclust','entrez','symbol']].to_csv('c1_network.csv')
gc2 = genedf[genedf['kclust']==2]
gc2[['tight10','kclust','entrez','symbol']].to_csv('c2_network.csv')

set1 = set(gc1['symbol'].tolist())
set2 = set(gc2['symbol'].tolist())
intersec = set1.intersection(set2)
int_list =['NFKB1',
'PRKCQ',
'RELA',
'RPS6KB1',
'IRS1',
'AR',
'HSPB1',
'JUN',
'CHEK1',
'IKBKB',
'RB1',
'IL6',
'IGF1R']
gc1_int = gc1[gc1['symbol'].isin(list(int_list))]
gc2_int = gc2[gc2['symbol'].isin(list(int_list))]
gg1 = gc1_int['symbol']
gg2 = gc2_int['symbol']

# 파일 및 변수 입력
MI_s = pd.read_csv('1. MI_s.csv',index_col=0)
MI_q = pd.read_csv('1. MI_q.csv',index_col=0)
MI_s = MI_s + pd.DataFrame.transpose(MI_s)
MI_q = MI_q + pd.DataFrame.transpose(MI_q)
MI_s_c1 = MI_s.loc[gc1_int.index,gc1_int.index]
MI_q_c2 = MI_q.loc[gc2_int.index,gc2_int.index]
MI_q_c1 = MI_q.loc[gc1_int.index,gc1_int.index]
MI_s_c2 = MI_s.loc[gc2_int.index,gc2_int.index]
