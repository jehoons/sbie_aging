import numpy as np
import numpy.linalg as la
import scipy as sp
import scipy.interpolate
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
        # regcoeff = np.polyfit(x, cond, 2)
        # reg = np.poly1d(regcoeff)
        # regdata = reg(t)
        tmp = scipy.interpolate.splrep(x, cond, k=1)
        regdata = scipy.interpolate.splev(t, tmp)
        intdata.append(regdata)
    intdatas.append(np.array(intdata))

# intdatas = conddatas#simulation with no interpolation

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
    T = Tp[:,0]
    for i in range(96):
        predT = R.dot(T) + T
        T=predT
    err = np.sqrt(np.mean((Tp[:,-1] - predT)**2))#RMSE
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


# import L100 list and make antibody list
import antibody_mapping as am

# os.chdir('L1000')
avg_up = pd.read_csv('L1000_average_up.csv', index_col = 0)
avg_down = pd.read_csv('L1000_average_down.csv', index_col = 0)
a375_up = pd.read_csv('L1000_A375_up.csv', index_col = 0)
a375_down = pd.read_csv('L1000_A375_down.csv', index_col = 0)
avg_up_list = avg_up['name'].values
avg_down_list = avg_down['name'].values
a375_up_list = a375_up['name'].values
a375_down_list = a375_down['name'].values
avg_up_ab = am.mapping(avg_up_list,ablist,genelist)
avg_down_ab = am.mapping(avg_down_list,ablist,genelist)
a375_up_ab = am.mapping(a375_up_list,ablist,genelist)
a375_down_ab = am.mapping(a375_down_list,ablist,genelist)
print(len(avg_up_ab))
print(len(avg_down_ab))
print(len(a375_up_ab))
print(len(a375_down_ab))
# make all combination
import allcombination as ac
avg_up_com = ac.c(avg_up_ab)
# print(avg_up_com)
avg_down_com = ac.c(avg_down_ab)
a375_up_com = ac.c(a375_up_ab)
a375_down_com = ac.c(a375_down_ab)
# print(avg_up_com[100])

#validating seneR using rapa data
#rapamycin is an inhibitor of mTOR which is known to inhibit senescence
file_name = 'avg_up_senR_result '
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(avg_up_com):
    print(u/len(avg_up_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = seneR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))
plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

file_name = 'avg_up_quiR_result'
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(avg_up_com):
    print(u/len(avg_up_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = quiR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))

plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)


###########222222222222222222222222222222222222222222222222222
file_name = 'avg_down_senR_result '
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(avg_down_com):
    print(u/len(avg_down_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = seneR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))
plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

file_name = 'avg_down_quiR_result'
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(avg_down_com):
    print(u/len(avg_down_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = quiR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))

plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

#####################################33333333333333333333333333333333
file_name = 'a375_up_senR_result '
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(a375_up_com):
    print(u/len(a375_up_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = seneR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))
plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

file_name = 'a375_up_quiR_result'
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(a375_up_com):
    print(u/len(a375_up_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = quiR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))

plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

###########444444444444444444444444444444444444
file_name = 'a375_down_senR_result '
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(a375_down_com):
    print(u/len(a375_down_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = seneR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))
plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

file_name = 'a375_down_quiR_result'
h = open('%s.txt'%file_name,'w')
s_err_total = []
q_err_total = []
for u,KD in enumerate(a375_down_com):
    print(u/len(a375_down_com))
    KDabidx = [ablist.index(ab) for ab in ablist if ab in KD]
    predrapaR = quiR.copy()
    for mab in KDabidx:
        predrapaR[:][mab] = 0
    # intrapa = intdatas[1]
    h.write('Knock down ab : ' )
    for item in KD:
        h.write("%s\t" % item)
    h.write('\n')
    h.write ('validation of modified senescence R to quiescence, and senescence data\n')
    # print('rappa: ',calcRerr_final_step(predrapaR, intrapa))
    q_err = calcRerr_final_step(predrapaR, intqui)
    s_err = calcRerr_final_step(predrapaR, intsene)
    h.write('pred vs qui : %s\n' %q_err)
    h.write('pred vs sen : %s\n' % s_err)
    q_err_total.append(q_err)
    s_err_total.append(s_err)
h.write('total qui err : %s\n' %(np.sum(q_err_total)))
h.write('total sen err : %s\n' %(np.sum(s_err_total)))

plt.figure()
plt.plot(s_err_total)
plt.savefig('%s_senErr.png'%file_name)
plt.figure()
plt.plot(q_err_total)
plt.savefig('%s_quiErr.png'%file_name)

    # # validating R with a random shuffled R matrix
    # valR = seneR.copy()
    # np.random.shuffle(valR)
    # print('validation using random shuffle of seneR')
    # # print('rappa: ',calcRerr_final_step(valR, intrapa))
    # print('qui :',calcRerr_final_step(valR, intqui))
    # print('sen :',calcRerr_final_step(valR, intsene))
# print('==================================')

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
plt.show()
'''
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

plt.show()
'''


