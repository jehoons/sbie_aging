import pandas as pd
import codecs

#import all the network files from previous analysis
#networkfiles = ['agingnetwork20p.sif', 'agingnetworkas13p.sif', 'agingnetworkia13p.sif', 'agingnetworkias16p.sif', 
#                'agingnetworkis14p.sif', 'agingnetworkqa13p.sif', 'agingnetworkqas16p.sif', 'agingnetworkqi13p.sif', 
#                'agingnetworkqia16p.sif', 'agingnetworkqias18p.sif', 'agingnetworkqis16p.sif', 'agingnetworkqr13p.sif', 
#                'agingnetworkqra16p.sif', 'agingnetworkqras18p.sif', 'agingnetworkqri16p.sif', 'agingnetworkqria18p.sif', 
#                'agingnetworkqris18p.sif', 'agingnetworkqrs15p.sif', 'agingnetworkqs12p.sif', 'agingnetworkra12p.sif', 
#                'agingnetworkras15p.sif', 'agingnetworkri13p.sif', 'agingnetworkria15p.sif', 'agingnetworkrias18p.sif', 
#                'agingnetworkris16p.sif', 'agingnetworkrs12p.sif']#network files for numerical cutoff of .5

#networkfiles = ['agingnetwork18p.sif', 'agingnetworkas12p.sif', 'agingnetworkia12p.sif', 'agingnetworkias15p.sif', 
#                'agingnetworkis12p.sif', 'agingnetworkqa12p.sif', 'agingnetworkqas14p.sif', 'agingnetworkqi11p.sif', 
#                'agingnetworkqia15p.sif', 'agingnetworkqias17p.sif', 'agingnetworkqis14p.sif', 'agingnetworkqr11p.sif', 
#                'agingnetworkqra14p.sif', 'agingnetworkqras16p.sif', 'agingnetworkqri14p.sif', 'agingnetworkqria16p.sif', 
#                'agingnetworkqris16p.sif', 'agingnetworkqrs14p.sif', 'agingnetworkqs10p.sif', 'agingnetworkra10p.sif', 
#                'agingnetworkras14p.sif', 'agingnetworkri11p.sif', 'agingnetworkria14p.sif', 'agingnetworkrias16p.sif', 
#                'agingnetworkris14p.sif', 'agingnetworkrs10p.sif']#network files for numerical cutoff of .75

networkfiles = ['agingnetwork13p.sif', 'agingnetworkas7p.sif', 'agingnetworkia7p.sif', 'agingnetworkias10p.sif', 
                'agingnetworkis7p.sif', 'agingnetworkqa7p.sif', 'agingnetworkqas10p.sif', 'agingnetworkqi6p.sif', 
                'agingnetworkqia10p.sif', 'agingnetworkqias12p.sif', 'agingnetworkqis9p.sif', 'agingnetworkqr7p.sif', 
                'agingnetworkqra10p.sif', 'agingnetworkqras12p.sif', 'agingnetworkqri10p.sif', 'agingnetworkqria12p.sif', 
                'agingnetworkqris11p.sif', 'agingnetworkqrs10p.sif', 'agingnetworkqs6p.sif', 'agingnetworkra6p.sif', 
                'agingnetworkras9p.sif', 'agingnetworkri6p.sif', 'agingnetworkria9p.sif', 'agingnetworkrias11p.sif', 
                'agingnetworkris9p.sif', 'agingnetworkrs6p.sif']#network files for numerical cutoff of .75

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
pd.DataFrame(netnetwork).to_csv('agingnetworkp_99.sif', sep='\t', index=None, header=None)




