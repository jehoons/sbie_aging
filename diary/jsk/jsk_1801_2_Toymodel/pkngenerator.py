import itertools as it

'''
This code generates numerous PKN_20 based on core PKN_19 and extravagant PKN_19 
'''
#import network information
corenetdata = []#list with core network information
corenetfile = open('PKN_19_core.sif', 'r')
for link in corenetfile:
    corenetdata.append(link.strip().split('\t'))
extranetdata = []#list with extravagant network information
extranetfile = open('PKN_19_extra.sif', 'r')
for link in extranetfile:
    extranetdata.append(link.strip().split('\t'))

#generate combination of variable links that are in the extravagant network but not in the core network
difflinks = [x for x in extranetdata if x not in corenetdata]
linkcomb = []#combination variable links from 0 links to all links
for n in range(len(difflinks) + 1):
    linkcomb = linkcomb + list(it.combinations(list(difflinks), n))

#generate list of all the possible pkn20s by concatenating variable links with the core network
pkn20s = []#list of all pkn20s; the links are in type sets because of itertools
for extralinks in linkcomb:
    pkn20s.append(corenetdata + list(extralinks))

#write each pkn_20 networks onto .sif files
for pknno, pkn20 in enumerate(pkn20s):
    pknno += 1
    filename = 'PKN_20_' + str(pknno) + '.sif'
    with open(filename, 'w') as f_out:#write file onto .sif file
        for link in pkn20:
            line = ''
            for item in link:
                line = line + item + '\t'
            line = line[:-1] + '\n'
            f_out.write(line)
