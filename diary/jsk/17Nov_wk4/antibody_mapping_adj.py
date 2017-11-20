#revising antibody mapping file to reflect the PKN_v14
import pandas as pd

#import relevant data
abmap = pd.read_csv('antibody_mapping_shell.csv', index_col=0)
agingnetwork = pd.read_csv('PKN_v15.csv')

#abmap shell contains the basic information of each antibodies
#and the mapped nodes; it only needs to fill the 'include network'
nabmap = abmap.copy()
for ab in abmap.index:
    if abmap['symbol'][ab] in list(agingnetwork['Source']) or abmap['symbol'][ab] in list(agingnetwork['Target']):
        nabmap['include network'][ab] = 'o'
    else:
        nabmap['include network'][ab] = 'x'

#export the resulting mapping file
nabmap.to_csv('antibody_mapping.csv')
