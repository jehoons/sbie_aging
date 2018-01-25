import numpy as np
import numpy.linalg as la
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

#import raw antibody data
ablist = []
with open('critablist.txt', 'r') as f:
    for ab in f:
        ablist.append(ab.strip('\n'))
data = pd.read_csv('rawdata.csv', index_col = 0)
data2 = data['Gene names']
data = data.drop(['Gene names'], axis=1)
#data = data.loc[ablist]
qui = np.array(data[['Quiescent 1h', 'Quiescent 6h', 'Quiescent 24h', 'Quiescent 96h']])
rapa = np.array(data[['Rapamycin 1h', 'Rapamycin 6h', 'Rapamycin 24h', 'Rapamycin 96h']])
igf = np.array(data[['IGF 0.5h', 'IGF 1h', 'IGF 3h', 'IGF 6h']])
apop = np.array(data[['Apoptosis 0.5h', 'Apoptosis 1h', 'Apoptosis 3h', 'Apoptosis 6h']])
sene = np.array(data[['Senescent 1h', 'Senescent 6h', 'Senescent 24h', 'Senescent 96h']])
conddatas = [qui, rapa, igf, apop, sene]
ablist = list(data.index)
genelist = list(data2.values)

#interpolate the given data
x1 = np.array([1, 6, 24, 96])#in hours; first type of time series for quiescent, rapamycin, and senescent
x2 = np.array([.5, 1, 3, 6])#in hours; second type of time series for igf and apoptosis
xs = [x1, x1, x2, x2, x1]
t1 = np.array(range(96 + 1))
t2 = np.array(range(6 * 2 + 1)) / 2
ts = [t1, t1, t2, t2, t1]
intdatas = []#[qui, rapa, igf, apop, sene]; interpolated dataset
for x, conddata, t in zip(xs, conddatas, ts):
    intdata = []
    for cond in conddata:
        regcoeff = np.polyfit(x, cond, 2)
        reg = np.poly1d(regcoeff)
        regdata = reg(t)
        #tmp = sp.interpolate.splrep(x, cond, k=1)
        #regdata = sp.interpolate.splev(t, tmp)
        intdata.append(regdata)
    intdatas.append(np.array(intdata))

intdatas = conddatas#simulation with no interpolation

#setting relevant constants
thr = 4#SVD threshold

def calcR(T, dT, thr):
    U, S, V = la.svd(T)
    tS = np.zeros((U.shape[0], V.shape[0]))#truncated S matrix
    
    # Abandon the eigen values whose arg is less than 'thr'.
    if thr >= S.size:#if threshold is set greater than the number of time points, it is fixed
        thr = S.size
    for i in range(thr):#write truncated S matrix
        tS[i, i] = S[i]
    
    # Get the inverse of tS with pseudo-inverse
    tS_inv = la.pinv(tS)
    
    # Infer matrix R
    R = np.dot(dT.dot(V.T), tS_inv.dot(U.T))
    return(R)

#Calculating R matrix for each conditions
Rs = []
for intdata in intdatas:
    T = intdata[:, :-1]
    dT = intdata[:, 1:] - intdata[:, :-1]
    R = calcR(T, dT, thr)
    Rs.append(R)
'''
#Calculating correlations between R matrices
R96 = [Rs[0].flatten(), Rs[1].flatten(), Rs[4].flatten()]
R6 = [Rs[2].flatten(), Rs[3].flatten()]
print(np.corrcoef(R96))
print(np.corrcoef(R6))
'''
#Calculating single composite R for quiescence and senescence
#       senescence  quiescence
#IGF1R  1           0
#AMPK   0           1
intqui = intdatas[0]
intsene = intdatas[4]
igf1rabs = ['IGF1R (Ab-1161)', 'IGF1R (Ab-1165/1166)', 'IGF1R (Phospho-Tyr1161)', 'IGF1R (Phospho-Tyr1165/1166)']
ampkabs = ['AMPK1/AMPK2 (Ab-485/491)', 'AMPK1 (Ab-172) ', 'AMPK beta1 (Ab-182)', 'AMPK1 (Phospho-Thr172)', 'AMPK1/AMPK2 (Phospho-Ser485/491)', 'AMPK beta1 (Phospho-Ser182)']
igf1rabidx = [ablist.index(ab) for ab in ablist if ab in igf1rabs]
ampkabidx = [ablist.index(ab) for ab in ablist if ab in ampkabs]
for qab in igf1rabidx:
    intqui[qab][:] = 0
for sab in ampkabidx:
    intsene[sab][:] = 0
qT = intqui[:, :-1]
sT = intsene[:, :-1]
dqT = intqui[:, 1:] - intqui[:, :-1]
dsT = intsene[:, 1:] - intsene[:, :-1]
T = np.concatenate((qT, sT), axis = 1)
dT = np.concatenate((dqT, dsT), axis = 1)
cR = calcR(T, dT, thr)#composite R for qui and sene
#R96.append(cR.flatten())
#print(np.corrcoef(R96))

#validating cR against quiR and seneR
quiR = Rs[0]#R for quiescence
seneR = Rs[4]#R for senescence
qscR = cR#R for both quiescence and senescence

def calcRerr(R, Tp):
    T = Tp[:, :-1]
    predT = R.dot(T) + T
    err = np.sqrt(np.mean((Tp[:, 1:] - predT)**2))#RMSE
    return err
def calcRerr_final_step(R, Tp):
    T = Tp[:, :-1]
    predT = R.dot(T) + T
    err = np.sqrt(np.mean((Tp[:, -1] - predT[:,-1])**2))#RMSE
    return err
#print('internal validation of quiescence and senescence R')
#print(calcRerr(quiR, intqui))
#print(calcRerr(seneR, intsene))
#print('validation of composite R on quiescence and senescence data')
#print(calcRerr(qscR, intqui))
#print(calcRerr(qscR, intsene))
#print('cross-error between qui and sene')
#print(calcRerr(quiR, intsene))
#print(calcRerr(seneR, intqui))

#validating seneR using rapa data
#rapamycin is an inhibitor of mTOR which is known to inhibit senescence

KD_list = ['SP3','AKR7A2','ATP1A3','BAK1','FZD5','ADAM15','SLC7A1','MAPK1','CORO1A','BRAF','COX4I1','DHPS','PPP2R5B','KRAS']
for gene in KD_list:
    intrapa = intdatas[1]
    KD=[]
    mtorabidx=[]
    x = [i for i,x in enumerate(genelist) if gene in str(x)]
    for i in x:
        KD.append(ablist[i])
    if not(KD):
        print('no match:',gene,KD)
        print('==================================')
        continue
    elif KD:
        print('match:',gene, KD)

    mtorabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = seneR.copy()
    for mab in mtorabidx:
        predrapaR[:][mab] = 0
    print('validation of predicted rapamycin R to rapamycin, quiescence, and senescence data')
    print(calcRerr_final_step(predrapaR, intrapa))
    print(calcRerr_final_step(predrapaR, intqui))
    print(calcRerr_final_step(predrapaR, intsene))
    # validating R with a random shuffled R matrix

    valR = seneR.copy()
    np.random.shuffle(valR)
    print('validation using random shuffle of seneR')
    print(calcRerr_final_step(valR, intsene))
    print(calcRerr_final_step(valR, intqui))
    print(calcRerr_final_step(valR, intrapa))
    print('==================================')
#need to reconstruct network by manually adjusting threshold hyperparameter
#Implementing Markov Cluster Algorithm on each R matrices
#Power of 2 and Inflation of 2
def normalize(M):#M is a two dimensional np.array() object
    return M / np.sum(M, axis = 0)
'''
#Markov Clustering Algorithm; an attempt to infer a network from the R matrix
def MCL(R, n):#R is the input matrix; two dimensional np.array() object; n is number of iterations
    if n == 0:
        print(R)
        return R
    else:
        print(R)
        R = R.dot(R)
        R = np.square(R)
        R = normalize(R)
        R[R < 1e-10] = 0
        return MCL(R, n - 1)

mclRs = []
for R in Rs:
    mclR = MCL(R, 6)
    mclRs.append(mclR)
'''
'''
#plotting to compare interpolated data to real data
idx_cond = 0
idx_plot = 2#97
plt.plot(range(97), intdatas[idx_cond][idx_plot], 'ro')
plt.plot(x1, qui[idx_plot])
plt.ylim(ymin=0)
'''

#associated network inference using sigma of R to the nth power
def sigRnth(R, n):
    sigR = np.zeros(np.shape(R))
    for i in range(1, n + 1):
        sigR = sigR + la.matrix_power(R, i)
    return sigR

tmpR = quiR
tmpR = sigRnth(tmpR, 8)
thres = .3#1.6
minset = [x for x in list(tmpR.flatten()) if x < -thres]
maxset = [x for x in list(tmpR.flatten()) if x > thres]
plt.hist(list(tmpR.flatten()), bins = np.arange(-thres, thres, .001), range = (-thres, thres))

# # plt.show()



