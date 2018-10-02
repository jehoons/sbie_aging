# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 11:07:09 2018

@author: ansu
"""

import pandas as pd
import scipy.stats
import numpy as np

#data_0 = pd.read_csv('0sample100.csv')
#data_1 = pd.read_csv('1sample100.csv')
#data_2 = pd.read_csv('2sample100.csv')
#
#
#
#data_0 = pd.read_csv('0sample500.csv')
#data_1 = pd.read_csv('1sample500.csv')
#data_2 = pd.read_csv('2sample500.csv')
df = pd.DataFrame()
for i in [100,200,300,400,500,600,700]:
    data_0 = pd.read_csv('0sample'+str(i)+'.csv')
    data_1 = pd.read_csv('1sample'+str(i)+'.csv')
    data_2 = pd.read_csv('2sample'+str(i)+'.csv')
    
    
    
    def naming(df):
        df["Unnamed: 1"].map(str)
        df = df.fillna(0)
        df["name"] = df["Unnamed: 0"] + df["Unnamed: 1"].map(str)
        return(df)
    data_0 = naming(data_0)        
    data_1 = naming(data_1)         
    data_2 = naming(data_2) 
    
    data_0 = data_0.sort_values(by=['count'], ascending=[False])
    data_0 = data_0.iloc[0:50,:]
    data_1 = data_1[data_1['name'].isin(data_0['name'])]
    data_2 = data_2[data_2['name'].isin(data_0['name'])]
    
    data_0 = data_0[data_0['name'].isin(data_1['name'])]
    data_2 = data_2[data_2['name'].isin(data_1['name'])]
    
    data_0 = data_0[data_0['name'].isin(data_2['name'])]
    data_1 = data_1[data_1['name'].isin(data_2['name'])]
    
    data_0 = data_0.sort_values(by=['name'], ascending=[False])
    data_1 = data_1.sort_values(by=['name'], ascending=[False])
    data_2 = data_2.sort_values(by=['name'], ascending=[False])
    
    
    
    c1 = data_0['count'].tolist()
    c2 = data_1['count'].tolist()
    c3 = data_2['count'].tolist()
    
    
    def calculate_rank(vector):
          a={}
          rank=1
          for num in sorted(vector, reverse=True):
              if num not in a:
                  a[num]=rank
                  rank=rank+1
          return[a[i] for i in vector]
    c1_rank = calculate_rank(c1)
    c2_rank = calculate_rank(c2)
    c3_rank = calculate_rank(c3)
    rho1,pval1 = scipy.stats.spearmanr(c1,c2)
    rho2,pval2 = scipy.stats.spearmanr(c1,c3)
    rho3,pval3 = scipy.stats.spearmanr(c2,c3)
    df[str(i)] = [rho1,rho2,rho3]
df.loc['average',:] = df.mean(axis=0)

import matplotlib.pyplot as plt
df_av = df.loc['average',:]
df_av.plot(kind='line')
plt.figure(); 
df_av.plot(kind='line')
plt.xlabel('The number of Model')
plt.ylabel('Average rank correlation coefficeint')
plt.title('Target rank correlation according to the number of model')
plt.ylim([0.8,1.0])
#plt.xlim([0,])
plt.show()
