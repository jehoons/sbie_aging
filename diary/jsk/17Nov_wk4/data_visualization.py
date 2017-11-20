#This code visualize the general trend of the data
#This is useful in determining general trend of the antibody data
#to select 'good' antibody for the midas data
import pandas as pd
import matplotlib.pyplot as plt

#import the relevant data
abmap = pd.read_csv('antibody_mapping.csv', index_col=0)
rawdata = pd.read_csv('raw_data.csv', index_col=0)

#exclude IGF 24h data for consensus
rawdata = rawdata.drop('IGF 24h', axis = 1)

#create dictionary of nodes and its respective antibodies
abmapdic = {}
for ab in abmap.index:
    if abmap['include network'][ab] == 'o':
        if abmap['symbol'][ab] in abmapdic:
            abmapdic[abmap['symbol'][ab]].append(ab)
        else:
            abmapdic[abmap['symbol'][ab]] = [ab]

#for each antibodies in abmapdic, convert raw data into MIDAS-like format
#and store it in midasdic; it has the same structure as abmapdic
midasdic = {}
for node in abmapdic.keys():
    midasdic[node] = []
    for ab in abmapdic[node]:
        abmidas = []
        for val in rawdata.loc[ab][1:]:
            abmidas.append(val / rawdata['control'][ab])
        midasdic[node].append([ab, abmidas])

#take a node name from the user to plot its antibodies
print(abmapdic.keys())
pltnode = input("Input node:\n")

#plot the antibody data for the given node; all of the conditions are separated individually
pltnum = 1
figsubplot = []
fig = plt.figure()
#show the plots individually
'''
colnum=5
rownum=len(midasdic[pltnode])
fig.set_figheight(2*rownum)
fig.set_figwidth(2*colnum)
for abdata in midasdic[pltnode]:
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][0:4])#plot for quiescence
    figsubplot[pltnum-1].set_title('Quiescence')
    figsubplot[pltnum-1].set_ylabel(abdata[0])
    pltnum += 1
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][4:8])#plot for mtor inhibition
    figsubplot[pltnum-1].set_title('mTOR Inhibition')
    pltnum += 1
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][8:12])#plot for IGF
    figsubplot[pltnum-1].set_title('IGF')
    pltnum += 1
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][12:16])#plot for apoptosis
    figsubplot[pltnum-1].set_title('Apoptosis')
    pltnum += 1
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][16:20])#plot for senescence
    figsubplot[pltnum-1].set_title('Senescence')
    pltnum += 1
'''
#show them all on one plot
colnum=1
rownum=len(midasdic[pltnode])
fig.set_figheight(4*rownum)
fig.set_figwidth(6*colnum)
for abdata in midasdic[pltnode]:
    figsubplot.append(fig.add_subplot(rownum, colnum, pltnum))
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][0:4], label='Quiescence')#plot for quiescence
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][4:8], label='mTOR Inhibition')#plot for mtor inhibition
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][8:12], label='IGF')#plot for IGF
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][12:16], label='apoptosis')#plot for apoptosis
    figsubplot[pltnum-1].plot([0,1,2,3,4], [1] + abdata[1][16:20], label='Senescence')#plot for senescence
    figsubplot[pltnum-1].set_ylabel(abdata[0])
    leg = figsubplot[pltnum-1].legend(loc='upper left')
    pltnum += 1
plt.tight_layout()
plt.show()














