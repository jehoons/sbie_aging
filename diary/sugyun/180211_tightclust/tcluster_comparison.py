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

# load data
MI_diff_complete = pd.read_csv('2. MI_diff_complete.csv',index_col=0)
Kclust = pd.read_csv('2. MI_diff_cluster.csv',index_col=0)
tightclust10 = pd.read_csv('tcluster1.txt',index_col=0,delimiter =' ')
tightclust5 = pd.read_csv('tcluster2.txt',index_col=0,delimiter =' ')
tightclust4 = pd.read_csv('tcluster3.txt',index_col=0,delimiter =' ')
# edit data
MI_diff_complete['tight10']=tightclust10[['x']].values
MI_diff_complete['tight5']=tightclust5[['x']].values
MI_diff_complete['tight4']=tightclust4[['x']].values
MI_diff_complete.sort_index(inplace=True)
Kclust.sort_index(inplace=True)
MI_diff_complete['kclust']=Kclust['clustering_label']
MI_diff_complete['entrez'] = Kclust['entrez']
result = pd.DataFrame.copy(MI_diff_complete)
result = result.sort_values(by=['kclust'],ascending=True)
result.to_csv('result.csv')

# save cluster and count the tight by kclust
clust_entrez = {}
for k in range(4):
    cnum = result[result['kclust']==k]
    cnum = cnum[cnum['tight10']>0]
    cnum = cnum.sort_values(by=['tight10'],ascending=True)
    gene_list = cnum['entrez'].values.flatten().tolist()
    gene_list_edited=[]
    for i in gene_list:
        if '/' in i:
            new = i.split('/')
            gene_list_edited.extend(new)
        else:
            new = i
            gene_list_edited.append(new)
    print(gene_list_edited)
    clust_entrez['c%s'%k] = set(gene_list_edited)
    cnum.to_csv('freq_tight10_kclust%s.csv' %k)
    freq = cnum['tight10'].value_counts()
    freq.sort_index(inplace=True)
    print(freq)
intersect = clust_entrez['c1'].intersection(clust_entrez['c2'])
for i in clust_entrez:
    with open('%s_tight.txt'%i,'w') as f:
        for j in clust_entrez[i]:
            if j in intersect:
                continue
            else:
                f.write('%s\n' %j)


# TSNE 그리기
#result = result[result['tight10']>0]
tight10 = result.pop('tight10')
tight5 = result.pop('tight5')
tight4 = result.pop('tight4')
kclust = result.pop('kclust')
entrez = result.pop('entrez')
#X = result.as_matrix()
#model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000,metric='euclidean')
##model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000)
#np.set_printoptions(suppress=True)
#X2 =model.fit_transform(X)
#x, y = X2.T
#plt.scatter(x, y, c=extra2.tolist(), cmap='jet')
##plt.legend(range(1,5),loc=4)
#plt.show()
#plt.scatter(x, y, c=clabel.tolist(), cmap='tab10')
##plt.legend(range(1,5),loc=4)
#plt.show()


# draw heatmap
#clabel = result.pop('clustering_label')
#entrez = MI_diff_complete.pop('entrez')
sns.set(color_codes=True)
lut = dict(zip(range(12), ['b','g','r','c','m','y','b','w','navy','aqua','gold','darkred']))


row_colors = tight10.map(lut)
g = sns.clustermap(result, cmap="mako", robust=True, row_colors = row_colors)