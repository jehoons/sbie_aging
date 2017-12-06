#This code converts raw data into MIDAS format
#It reads raw_data and antibody_mapping and selects
#the mapped antibodies and reformat its data into MIDAS
import pandas as pd
from collections import OrderedDict
import codecs
from itertools import product

#import the relevant data
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)
rawdata = pd.read_csv('raw_data.csv', index_col=0)
midas = pd.read_csv('MIDAS_v16_shell.csv', index_col=None)

#exclude IGF 24h data for consensus
rawdata = rawdata[['control', 'Quiescent 96h', 'IGF 24h', 'Apoptosis 6h', 'Senescent 96h']]

#create dictionary of nodes and its respective antibodies
abmapdic = OrderedDict()
for ab in abmap.index:
    if abmap['include network'][ab] == 'o':
        if abmap['symbol'][ab] in abmapdic:
            abmapdic[abmap['symbol'][ab]].append(ab)
        else:
            abmapdic[abmap['symbol'][ab]] = [ab]
abmapdic['TP53_p'] = abmapdic['TP53']#minor accomodation for TP53_p; temporary measure

#set time points on midas data for nodes represented by antibody data
for node in abmapdic.keys():
    abcond = []
    abcond.append(0)
    for i in range(4):
        abcond.append(1)
    colname = 'DA:' + node
    midas[colname] = abcond

#set empty columns representing values for each nodes
for node in abmapdic.keys():
    colname = 'DV:' + node
    midas[colname] = [0] * len(midas.index)

#based on the antibodies found, the relevant antibodies were chosen manually
#the all combination of antibody combination possible are converted to MIDAS format
#manselab=pd.read_csv('manual_antibody_selection.csv', index_col=0)
manselab = []
with codecs.open('manual_antibody_selection.txt', "r", encoding="utf-8") as f_in:
    for line in f_in:
        item = line.strip('\n').strip('\r').split('\t')
        manselab.append(item)
mselabdic = OrderedDict()
for node in manselab:
    if len(node) == 1:
        pass
    else:
        mselabdic[node[0]] = node[1:]
mselabdic['TP53_p'] = mselabdic['TP53']#minor accomodation for TP53_p; temporary measure

#for each antibodies in mselabdic, convert raw data into MIDAS format
#and store it in midasdic; it has the same structure as mselabdic
midasdic = OrderedDict()
for node in mselabdic.keys():
    midasdic[node] = []
    for ab in mselabdic[node]:
        abmidas = []
        abmidas.append(rawdata['control'][ab])
        for val in rawdata.loc[ab][1:]:
            abmidas.append(val)
        midasdic[node].append([ab, abmidas])

#make a preprocessed boolean form of midasdic;
#The control are considered to be at 0
#if the data value is less than control -> 0
#if the data value is greater than control -> 1
midasdicbool = midasdic
for node in midasdicbool.keys():
    for abmidas in midasdicbool[node]:
        for idx, val in enumerate(abmidas[1][1:]):
            if val > abmidas[1][0]:
                abmidas[1][idx + 1] = 1
            else:
                abmidas[1][idx + 1] = 0
        abmidas[1][0] = 1

'''
#unzip midasdic in order to make individual MIDAS format data
#export the MIDAS file
midasset = list(product(*midasdic.values()))
abused = []
for n, tmidas in enumerate(midasset):
    tabused = []
    for nodename, node in zip(midasdic.keys(), tmidas):
        colname = 'DV:' + nodename
        midas[colname] = node[1]
        tabused.append(node[0])
    midas.to_csv('MIDAS_v16_' + str(n) + '.csv', index=False)
    abused.append(tabused)
abused = pd.DataFrame(abused, columns = midasdic.keys())
abused.to_csv('abused.csv')
'''













