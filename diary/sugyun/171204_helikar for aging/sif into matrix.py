import pandas as pd
import numpy as np
f = open('PKN_v17.txt','r')
g = open('PKN_v17_weighted.txt','w')
node_list = []
for line in f:
    line = line.strip().split('\t')
    node_list.append(line[0])
    node_list.append(line[2])
node_set = set(node_list)
node_set = sorted(node_set)

matrix = np.zeros((len(node_set),len(node_set)))
f = open('PKN_v17.txt','r')
for line in f:
    line = line.strip().split('\t')
    row = node_set.index(line[0])
    col = node_set.index(line[2])
    matrix[col][row] = round(int(line[1]))
# print(matrix)
    # g.write("%s\t%s\n"%(line[0],line[1]))

np.savetxt('PKN_v17_weighted.txt', matrix.astype(int), fmt ='%.0f\t')

print(node_set)
print(len(node_set))