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
import mygene
import networkx as nx


# 파일 및 변수 입력
entrez = pd.read_csv('2. entrez.csv',index_col=0)
q = pd.read_csv('1. MI_q.csv',index_col=0)
s = pd.read_csv('1. MI_s.csv',index_col=0)
p = pd.read_csv('1. MI_p.csv',index_col=0)

# default 생성
mg = mygene.MyGeneInfo()
G = nx.Graph()
col_num_q = q.shape[0]

### symbol로 치환하기
#entrez['symbol'] = entrez['GeneID']
#for n, ent in enumerate(entrez['GeneID'].values.tolist()):
#    if "/" in ent:
#        a = ent.split('/')
#    else:
#        a = [ent]
#    a_g = []
#    for a_ent in a:
#        try:
#            a_g.append(mg.getgene(a_ent,'symbol')['symbol'])
#        except:
#            continue
#    sym = '/'.join( a_g)
#    entrez.iloc[n,1]=sym
#entrez.to_csv('entrezwithsym.csv')

# cutoff graph 만들기 
entrez = pd.read_csv('entrezwithsym.csv',index_col=0)
entrez = entrez['symbol']
for n,i in enumerate(range(col_num_q)):
    q_c = pd.DataFrame(q.ix[i,:])
    c = q_c[q_c.gt(1.2)]
    r = c.columns.tolist()
    c = c.dropna(subset = [r])
    a = entrez.ix[n,'symbol']
    if "/" in a:
        a = a.split('/')
    else:
        a = [a]
    if a == []:
        continue
    for j in c.index.tolist():
        b = entrez.ix[j,'symbol']
        if "/" in b:
            b = b.split('/')
        else:
            b=[b]
        if b == []:
            continue
        for a_sym in a:
            for b_sym in b:
                G.add_edge(a_sym,b_sym)


# omnipath 만들기
H = nx.Graph()
f = open('180201_omnipath_symbol.txt','r')
for line in f:
    line = line.strip().split('\t')
    H.add_edge(line[0],line[2])

## intersection edges (the node must be same) 
R=G.copy()
R.remove_nodes_from(n for n in G if n not in H)
T = H.copy()
T.remove_nodes_from(n for n in H if n not in G)
Q = nx.intersection(R,T)


###
q = s
mg = mygene.MyGeneInfo()
G = nx.Graph()
col_num_q = q.shape[0] 

entrez = pd.read_csv('entrezwithsym.csv',index_col=0)
entrez = entrez['symbol']
for n,i in enumerate(range(col_num_q)):
    q_c = pd.DataFrame(q.ix[i,:])
    c = q_c[q_c.gt(1.2)]
    r = c.columns.tolist()
    c = c.dropna(subset = [r])
    a = entrez.ix[n,'symbol']
    if "/" in a:
        a = a.split('/')
    else:
        a = [a]
    if a == []:
        continue
    for j in c.index.tolist():
        b = entrez.ix[j,'symbol']
        if "/" in b:
            b = b.split('/')
        else:
            b=[b]
        if b == []:
            continue
        for a_sym in a:
            for b_sym in b:
                G.add_edge(a_sym,b_sym)


# omnipath 만들기
H = nx.Graph()
f = open('180201_omnipath_symbol.txt','r')
for line in f:
    line = line.strip().split('\t')
    H.add_edge(line[0],line[2])

## intersection edges (the node must be same) 
R=G.copy()
R.remove_nodes_from(n for n in G if n not in H)
T = H.copy()
T.remove_nodes_from(n for n in H if n not in G)
S = nx.intersection(R,T)

###
q = p
mg = mygene.MyGeneInfo()
G = nx.Graph()
col_num_q = q.shape[0] 


entrez = pd.read_csv('entrezwithsym.csv',index_col=0)
entrez = entrez['symbol']
for n,i in enumerate(range(col_num_q)):
    q_c = pd.DataFrame(q.ix[i,:])
    c = q_c[q_c.gt(1.2)]
    r = c.columns.tolist()
    c = c.dropna(subset = [r])
    a = entrez.ix[n,'symbol']
    if "/" in a:
        a = a.split('/')
    else:
        a = [a]
    if a == []:
        continue
    for j in c.index.tolist():
        b = entrez.ix[j,'symbol']
        if "/" in b:
            b = b.split('/')
        else:
            b=[b]
        if b == []:
            continue
        for a_sym in a:
            for b_sym in b:
                G.add_edge(a_sym,b_sym)


# omnipath 만들기
H = nx.Graph()
f = open('180201_omnipath_symbol.txt','r')
for line in f:
    line = line.strip().split('\t')
    H.add_edge(line[0],line[2])

## intersection edges (the node must be same) 
R=G.copy()
R.remove_nodes_from(n for n in G if n not in H)
T = H.copy()
T.remove_nodes_from(n for n in H if n not in G)
P = nx.intersection(R,T)

## save edge
def edge_save(L, filename):
    f = open(filename,'w')
    for i,o in L.edges:
        f.write('%s\t%s\n' %(i,o))
    f.close()
edge_save(S,'s.txt')
edge_save(Q,'q.txt')
edge_save(P,'p.txt')


P2=P.copy()
S2=S.copy()
Q2=Q.copy()
P2.remove_nodes_from(n for n in P if n not in S)
S2.remove_nodes_from(n for n in S if n not in P)
PS = nx.intersection(P2,S2)
edge_save(PS,'PS.txt')

P2=P.copy()
S2=S.copy()
Q2=Q.copy()
Q2.remove_nodes_from(n for n in Q if n not in S)
S2.remove_nodes_from(n for n in S if n not in Q)
SQ = nx.intersection(Q2,S2)
edge_save(SQ,'SQ.txt')

P2=P.copy()
S2=S.copy()
Q2=Q.copy()
P2.remove_nodes_from(n for n in P if n not in Q)
Q2.remove_nodes_from(n for n in Q if n not in P)
PQ = nx.intersection(P2,Q2)
edge_save(PQ,'PQ.txt')


P3=PQ.copy()
S3 = SQ.copy()
P3.remove_nodes_from(n for n in PQ if n not in SQ)
S3.remove_nodes_from(n for n in SQ if n not in PQ)
SPQ = nx.intersection(P3,S3)
edge_save(SPQ,'SPQ.txt')





