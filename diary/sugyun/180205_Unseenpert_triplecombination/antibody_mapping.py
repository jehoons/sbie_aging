import numpy as np
import numpy.linalg as la
import scipy as sp
import scipy.interpolate
import pandas as pd
import matplotlib.pyplot as plt

def mapping(List,ablist,genelist):
    ab = []
    for gene in List:
        try:
            x = [i for i,x in enumerate(genelist) if gene in str(x)]
        except:
            print(gene)
        for i in x:
            ab.append(ablist[i])
    return(ab)