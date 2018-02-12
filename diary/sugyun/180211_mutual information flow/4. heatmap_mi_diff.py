import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from sklearn.manifold import TSNE
from matplotlib import cm

# load data
MI_diff_complete = pd.read_csv('2. MI_diff_cluster.csv',index_col=0)

## save cluster and pickle
#clust_entrez = {}
#for k in range(4):
#    cnum = MI_diff_complete[MI_diff_complete['clustering_label']==k]
#    gene_list = cnum['entrez'].values.flatten().tolist()
#    gene_list_edited=[]
#    for i in gene_list:
#        if '/' in i:
#            new = i.split('/')
#            gene_list_edited.extend(new)
#        else:
#            new = i
#            gene_list_edited.append(new)
#    clust_entrez['cluster'+str(k)] = gene_list_edited
#with open('4. pickle.pickle', 'wb') as f:
#        pickle.dump(clust_entrez, f)
# load cluster pickle
with open('4. pickle.pickle','rb') as f:
    dat = pickle.load(f)
    for key in dat.keys():
        exec(key+'= dat[key]')
with open('4. cluster1_high in sen.txt','w') as g:
    for i in cluster1:
        g.write('%s\n'%i)
with open('4. cluster2_low in sen.txt','w') as g:
    for i in cluster2:
        g.write('%s\n'%i)
    
    

## draw heatmap
#clabel = MI_diff_complete.pop('clustering_label')
#entrez = MI_diff_complete.pop('entrez')
#sns.set(color_codes=True)
#lut = dict(zip(range(4), ['navy','aqua','gold','darkred']))
#row_colors = clabel.map(lut)
#g = sns.clustermap(MI_diff_complete, cmap="mako", robust=True, row_colors = row_colors)

## TSNE 그리기
#X = MI_diff_complete.as_matrix()
#model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=1000,metric='euclidean')
#np.set_printoptions(suppress=True)
#X2 =model.fit_transform(X)
#x, y = X2.T
#plt.scatter(x, y, c=clabel.tolist(), cmap='jet')
#plt.show()