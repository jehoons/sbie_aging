# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 16:58:17 2018

@author: ansu
"""

import seaborn as sns; sns.set(color_codes=True)
iris = sns.load_dataset("iris")
species = iris.pop("species")
g = sns.clustermap(iris)