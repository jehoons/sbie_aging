import pandas as pd
import codecs

#import all the network files from previous analysis
networkfiles = ['agingnetwork20p.sif', 'agingnetworkas13p.sif', 'agingnetworkia13p.sif', 'agingnetworkias16p.sif', 
                'agingnetworkis14p.sif', 'agingnetworkqa13p.sif', 'agingnetworkqas16p.sif', 'agingnetworkqi13p.sif', 
                'agingnetworkqia16p.sif', 'agingnetworkqias18p.sif', 'agingnetworkqis16p.sif', 'agingnetworkqr13p.sif', 
                'agingnetworkqra16p.sif', 'agingnetworkqras18p.sif', 'agingnetworkqri16p.sif', 'agingnetworkqria18p.sif', 
                'agingnetworkqris18p.sif', 'agingnetworkqrs15p.sif', 'agingnetworkqs12p.sif', 'agingnetworkra12p.sif', 
                'agingnetworkras15p.sif', 'agingnetworkri13p.sif', 'agingnetworkria15p.sif', 'agingnetworkrias18p.sif', 
                'agingnetworkris16p.sif', 'agingnetworkrs12p.sif']

agnetlist = []#contains all of the calculated network
for netfile in networkfiles:
    agnet = []
    with codecs.open(netfile, "r", encoding="utf-8") as f_in:
        for line in f_in:
            interaction = line.strip().split()
            agnet.append(interaction)
    agnetlist.append(agnet)

#combine all of the networks into one
netnetwork = []
for agnet in agnetlist:
    for interaction in agnet:
        if not interaction in netnetwork:
            netnetwork.append(interaction)

#export the net network files to .sif format
pd.DataFrame(netnetwork).to_csv('agingnetworkp.sif', sep='\t', index=None, header=None)




