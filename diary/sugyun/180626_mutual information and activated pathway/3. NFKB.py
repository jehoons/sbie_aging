import pandas as pd

"""
nfkb의 세 pathway가 잘 클러스터링 되는지 확인한 결과
0번 클러스터에는 p105가
1번 클러스터에는 p65가 위치하고 있음을 확인함
(p100은 어디에 특별히 클러스터링 되지 않음)

"""
data = pd.read_csv('2. MI_diff_cluster.csv',index_col=0)
clus = data['clustering_label']

bool_list = []
for i,ab in enumerate(data.index):
    if 'NFkB' in ab:
        bool_list.append(i)
nfkb = clus.loc[data.index[bool_list]]
nfkbc0 = nfkb[nfkb==0]
nfkbc1 = nfkb[nfkb==1]
nfkbc2 = nfkb[nfkb==2]
nfkbc3 = nfkb[nfkb==3]
###############################################################################
"""
그렇다면 0번 클러스터와 1번 클러스터에 속한 NFKB의 각 sen과 qui에서의
그래프 양상을 보기로 함
"""
# Make a data frame
def makeDF(my_dataframe, filename):
    dic = {'x':range(4)}
    for i in my_dataframe.index:
        dic[i]=my_dataframe.loc[i,]   
    df=pd.DataFrame(dic)
    # style
    plt.style.use('seaborn-darkgrid')
     
    # create a color palette
    palette = plt.get_cmap('Set1')
     
    # multiple line plot
    num=0
    for column in df.drop('x', axis=1):
        num+=1
        plt.plot(df['x'], df[column], marker='', color=palette(num), linewidth=1, alpha=0.9, label=column)
     
    # Add legend
#    plt.legend(loc=1, ncol=2)
     
    # Add titles
    plt.title(filename, loc='left', fontsize=12, fontweight=0, color='orange')
    plt.xlabel("Time")
    plt.ylabel("fold change")
    plt.savefig('3. %s.png' % filename)
    plt.clf()
    return
    
data2 = pd.read_csv('1. fold_change.csv',index_col=0)
data2 = data2.iloc[:,1:]
data_sq = data2.iloc[:,[0,1,2,3,17,18,19,20]]
data_sq_nfkb = data_sq.loc[data.index[bool_list]]
c0 = data_sq_nfkb.loc[nfkbc0.index]
c1 = data_sq_nfkb.loc[nfkbc1.index]
c2 = data_sq_nfkb.loc[nfkbc2.index]
c3 = data_sq_nfkb.loc[nfkbc3.index]


# libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
 
# data
c0_q=c0.iloc[:,[0,1,2,3]]
c0_s=c0.iloc[:,[4,5,6,7]]
c1_q=c1.iloc[:,[0,1,2,3]]
c1_s=c1.iloc[:,[4,5,6,7]]
c2_q=c2.iloc[:,[0,1,2,3]]
c2_s=c2.iloc[:,[4,5,6,7]]
c3_q=c3.iloc[:,[0,1,2,3]]
c3_s=c3.iloc[:,[4,5,6,7]]
a = makeDF(c0_q,'co_q')
a = makeDF(c0_s,'c0_s')
a = makeDF(c1_q,'c1_q')
a = makeDF(c1_s,'c1_s')
a = makeDF(c2_q,'c2_q')
a = makeDF(c2_s,'c2_s')
a = makeDF(c3_q,'c3_q')
a = makeDF(c3_s,'c3_s')