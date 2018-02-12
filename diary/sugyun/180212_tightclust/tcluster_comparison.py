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
MI_s = pd.read_csv('MI_s_clustering.csv',index_col=0)
MI_q = pd.read_csv('MI_q_clustering.csv',index_col=0)

# data 수정
extra = result.pop('tight5')
extra = result.pop('tight4')
result.sort_index(inplace=True)
MI_s.sort_index(inplace=True)
MI_q.sort_index(inplace=True)
result['sclust'] = MI_s['clustering_label']
result['qclust'] = MI_q['clustering_label']
result = result.sort_values(by=['kclust'],ascending=True)

#각 클러스터의 각 tight clust 저장
Cresult = result[result['tight10']>0]
for i in range(4):
    kc = Cresult[Cresult['kclust']==i]
    for j in range(1,11):
        tc = kc[kc['tight10']==j]
        if tc.shape[0] == 0:
            continue
        tc[['kclust','tight10','entrez','symbol']].to_csv('C%s_T%s.csv'%(i,j))


# TSNE 그리기
result = result[result['tight10']>0]
tight10 = result.pop('tight10')
kclust = result.pop('kclust')
entrez = result.pop('entrez')
sym = result.pop('symbol')
slabel = result.pop('sclust')
qlabel = result.pop('qclust')


X = result.as_matrix()
model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000,metric='euclidean')
#model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000)
np.set_printoptions(suppress=True)
X2 =model.fit_transform(X)
x, y = X2.T
plt.scatter(x, y, c=kclust.tolist(), cmap='jet')
plt.show()
plt.scatter(x, y, c=tight10.tolist(), cmap='tab10')
plt.show()


#slabel2 = MI_s.pop('clustering_label')
#X = MI_s.as_matrix()
#model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000,metric='euclidean')
##model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000)
#np.set_printoptions(suppress=True)
#X2 =model.fit_transform(X)
#x, y = X2.T
#plt.scatter(x, y, c=slabel2.tolist(), cmap='jet')
#plt.show()

#qlabel2 = MI_q.pop('clustering_label')
#X = MI_q.as_matrix()
#model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000,metric='euclidean')
##model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000)
#np.set_printoptions(suppress=True)
#X2 =model.fit_transform(X)
#x, y = X2.T
#plt.scatter(x, y, c=qlabel2.tolist(), cmap='jet')
#plt.show()
