import json
from os.path import exists
from boolean3_addon import attr_cy
from sbie_aging.results import table_s1
from os.path import dirname
from ipdb import set_trace
import pandas as pd 
import itertools 

file_c = dirname(table_s1.__file__) + '/c/attractors.json'

file_d1 = dirname(table_s1.__file__) + '/d/point.csv'
file_d2 = dirname(table_s1.__file__) + '/d/point_short.csv'

with open(file_c, 'r') as fobj_c: 
    data = json.load(fobj_c)

dfpoint = pd.DataFrame([], columns=data['labels']+['key','ratio'])
k = 0 
for attr in data['attractors']:
    thisattr = data['attractors'][attr]
    if thisattr['type'] == 'point': 
        value = thisattr['value']
        bitvec = data['state_key'][value]
        # print(bitvec)
        # set_trace(); break
        dictdata = dict([(a,b) for a,b in zip(data['labels'], bitvec)])
        dfpoint.loc[k, data['labels']] = dictdata
        dfpoint.loc[k, 'key'] = value
        dfpoint.loc[k, 'ratio'] = thisattr['ratio']
        k += 1

heteroprofile = [] 
for label in data['labels']: 
    n = len(dfpoint.groupby(label).groups.keys())
    if n > 1 : 
        heteroprofile.append(label)

dfpointshort = dfpoint[heteroprofile+['key','ratio']]

dfpoint.to_csv(file_d1)
dfpointshort.to_csv(file_d2)
