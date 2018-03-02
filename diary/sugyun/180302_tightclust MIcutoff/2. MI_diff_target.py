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
MI_s = pd.read_csv('1. MI_s.csv',index_col=0)
MI_q = pd.read_csv('1. MI_q.csv',index_col=0)

# clustering
MI_s_complete = MI_s.add(pd.DataFrame.transpose(MI_s), fill_value=0)
MI_q_complete = MI_q.add(pd.DataFrame.transpose(MI_q), fill_value=0)

mat = MI_s_complete.as_matrix()
km = KMeans(n_clusters=4)
km.fit(mat)
# Get cluster assignment labels
labels = km.labels_

# Format results as a DataFrame
MI_s_complete['clustering_label'] = labels
MI_s_complete.to_csv('MI_s_clustering.csv')

mat = MI_q_complete.as_matrix()
km = KMeans(n_clusters=4)
km.fit(mat)
# Get cluster assignment labels
labels = km.labels_

# Format results as a DataFrame
MI_q_complete['clustering_label'] = labels
MI_q_complete.to_csv('MI_q_clustering.csv')