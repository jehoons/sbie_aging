import pandas as pd

data = pd.read_csv('gse2487_FCupdown.csv', index_col=0)
f = open('upatsen_list.txt', 'w')
data_up = data[data['logFC']>0]
data_down = data[data['logFC']<0]
up_list = data_up['Gene.symbol']
down_list = data_down['Gene.symbol']
up_set = set(up_list.values)
down_set = set(down_list.values)
print(up_set)
g = open('up_genesymbol.txt','w')
h = open('down_genesymbol.txt','w')
for i in up_set:
    if str(i) == 'nan':
        continue
    elif "///" in str(i):
        a = i.split('///')
        for k in a:
            g.write('%s\n'%(k))
    else:
        g.write('%s\n'%(i))
for i in down_set:
    if str(i) == 'nan':
        continue
    elif "///" in str(i):
        a = i.split('///')
        for k in a:
            h.write('%s\n'%(k))
    else:
        h.write('%s\n'%(i))
