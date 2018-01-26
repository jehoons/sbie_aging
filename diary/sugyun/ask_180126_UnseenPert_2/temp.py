import pandas as pd

data = pd.read_csv('rawdatafc.csv', index_col = 0)
nodelist = []
with open('nodelist.txt', 'r') as nodelistf:
    for node in nodelistf:
        nodelist.append(node.strip())

critab = []
ablist = [ab for ab in data.index]
for node in nodelist:
    for ab in ablist:
        if node in ab:
            critab.append(ab)

with open('critablist.txt', 'w') as f:
    for ab in critab:
        f.write(ab + '\n')
        