import igraph 
from igraph import *
import pandas as pd 
# import matplotlib.pyplot as plt
networkfile = 'network_2.txt'

df0 = pd.read_csv(networkfile)
g = Graph()

g.add_vertices(df0['#node1'].values.tolist())
g.add_vertices(df0['node2'].values.tolist())

for row in df0.iterrows(): 
    # print(row)
    node1 = row[1]['#node1']
    node2 = row[1]['node2']
    score = row[1]['combined_score']
    g.add_edge(node1, node2, score=score)


xxx
# nodepairs = df0[['#node1','node2']].values.tolist() 

# nodes = [] 
# for nodepair in nodepairs:
#     nodes = nodes + nodepair 

# nodes = [x for x in set(nodes)]
# node_map =

# for i,n in enumerate(nodes):



# xxx

# g = Graph(nodelist)

# xxx

# g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
# g.vs["name"] = ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"]
# g.vs["age"] = [25, 31, 18, 47, 22, 23, 50]
# g.vs["gender"] = ["f", "m", "f", "m", "f", "m", "m"]
# g.es["is_formal"] = [False, False, True, True, True, False, True, False, False]
# g["date"] = "2009-01-10"

# g.edge_betweenness()

# layout = g.layout("kk")

# plot(g,"output.png", layout=layout)

