### https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from sklearn.metrics import mutual_info_score
import numpy as np
from itertools import combinations

# 파일 및 변수 입력
entrez = pd.read_csv('2. entrez.csv',index_col=0)
q = pd.read_csv('1. MI_a.csv',index_col=0)
s = pd.read_csv('1. MI_s.csv',index_col=0)
col_count = q.shape[1]

# MI의 차이 구하기
diff = pd.DataFrame(pd.DataFrame(np.zeros((col_count,col_count)),columns=list(q.index),index=list(q.index)))
diff = s-q
#diff = abs(q-s)
diff.to_csv("2. MI_diff.csv")

# clustering
diff = pd.read_csv('2. MI_diff.csv',index_col=0)
diff_complete = diff.add(pd.DataFrame.transpose(diff), fill_value=0)
diff_complete.to_csv("2. MI_diff_complete.csv")
mat = diff_complete.as_matrix()
km = KMeans(n_clusters=4)
km.fit(mat)
# Get cluster assignment labels
labels = km.labels_

# Format results as a DataFrame
diff_complete['clustering_label'] = labels
diff_complete['entrez'] = entrez.values
diff_complete=diff_complete.sort_values(by=['clustering_label'],ascending=True)
diff_complete.to_csv("2. MI_diff_cluster.csv")

for i in range(4):
    diff_cluster = diff_complete[diff_complete['clustering_label']==i]
    gene_list = diff_cluster['entrez'].values.flatten().tolist()
    gene_list_edited=[]
    for i in gene_list:
        if '/' in i:
            new = i.split('/')
            gene_list_edited.extend(new)
        else:
            new = i
            gene_list_edited.append(new)
    print(gene_list_edited)

# TSNE 그리기
X = mat
model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=5000,metric='euclidean')
np.set_printoptions(suppress=True)
X2 =model.fit_transform(X)
x, y = X2.T
# colors = cm.rainbow(np.linspace(0, 1, 10))
plt.scatter(x, y, c=labels, cmap=cm.rainbow)
plt.savefig('2. MI_cluster.png')
plt.show()

"""
#################################################raw_data clustering and visualize######################################
fc = pd.read_csv('fold_change.csv',index_col=0)
### clustering
mat = fc.as_matrix()
km = KMeans(n_clusters=6)
km.fit(mat)
# Get cluster assignment labels
labels = km.labels_
# Format results as a DataFrame
fc['clustering_label'] = labels
fc['entrez'] = entrez.values
fc=fc.sort_values(by=['clustering_label'],ascending=False)
fc.to_csv("MI_fc_cluster.csv")
for i in range(6):
    new = fc[fc['clustering_label']==(i+1)]
    gene_list = new['entrez'].values.flatten().tolist()
    gene_list_edited=[]
    for i in gene_list:
        if '/' in i:
            new = i.split('/')
            gene_list_edited.extend(new)
        else:
            new = i
            gene_list_edited.append(new)
    print(gene_list_edited)

X = mat
model = TSNE(n_components=2, random_state=0,verbose=2,perplexity=100, n_iter=5000,metric='euclidean')
np.set_printoptions(suppress=True)
X2 =model.fit_transform(X)
# print(X2)
x, y = X2.T

# colors = cm.rainbow(np.linspace(0, 1, len(C)))
plt.scatter(x, y, c=labels, cmap=cm.rainbow)
plt.show()
"""